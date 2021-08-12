/*
 *
 * Dregion class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * Creates a hexgrid, divides it into Dirichlet domains
 * then solves KS equations and can morph to curved
 * boundaries
 *
 */
#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/ReadCurves.h>
#include <morph/HdfData.h>
#include <morph/BezCurve.h>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <hdf5.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifndef HEX_GEOM
#define HEX_GEOM
#include "hexGeometry.h"
#endif
#include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::ReadCurves;
using morph::Tools;
using morph::BezCurve;
using morph::BezCurvePath;
using morph::BezCoord;
using namespace std;

#ifndef NO_N
#define NO_N
#define nE(hi) (this->Hgrid->d_ne[hi])
#define HAS_nE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define nW(hi) (this->Hgrid->d_nw[hi])
#define HAS_nW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define nNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_nNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define nNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_nNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define nSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_nSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define nSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_nSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)
#endif

/*!
 * This class is used to build a HexGrid and then to dissect it
 * into Voroni regions based on a set of points centres currently
 * incorporated as a header file. The class contains
 * methods for identifying and dissecting the boundary into edges
 * and creating structures to hold vertices and edges and to do
 * the bookeeping necessary to keep track of adjacency. It also
 * contains methods to sectorize regions by angular or radial
 * sections. It contains methods to derive the Pearson r coefficient
 * for the vectors of values on edges. There are separate routines
 * for adjacent or randomly chosen edges. It contains Runge Kutta
 * solvers for the Keller-Segel equations, though these should be
 * delegated to the ksSolver class. It contains routines to derive
 * the area and perimeter of Voroni regions. It contains methods
 * to round the corners of the original Voroni tessellation and also
 * subsequent morphed tessellations.
 */

class DRegion
{
public:
    /*!
     * list of scalar objects at public scope
     */
    int scale; //scale of the initial HexGrid
    FLT xspan; //width of the HexGrid
    int n; //size of the initial HexGrid/
    FLT ds; //hex to hex distance
    string logpath; //where the analysis and restart files live
    const int base = 1000; //for hashing of (i,j) pairs for each edge
    const FLT tan30 = 0.5773502692;
    FLT hexArea;
    unsigned int NUMPOINTS;
    /*
     * list of containers at public scope
     */
    vector<vector<int> > N; // hex neighbourhood
    vector<morph::Hex> vHexInGrid; // hexen for the original HexGrid as a vector rather than a list
    vector<vector<int> > region; //for each hex sorted list of regions as ints
    vector<bool> innerRegion;
    vector<list<morph::Hex>> regionHex; //for each region list of hexes it contains
    vector<vector<int>> hexRegionList; //for each hex neighbour regions as ints
    vector<vector<int>> regionList; //for each region, regions that are its neighbours
    vector<vector<int>> regionVertex; //for each region, vertices that bound it sorted by angle
    vector<vector<int>> sortedBoundary; // for each reagion indices of boundary or cutter hexes sorted by angle
    vector<vector<FLT>> sortedBoundaryPhi; // for each region phi values hexes sorted by angle
    vector<vector<FLT>> sortedBoundaryNN; // for each region NN values sorted by angle
    vector<int> edgeIndex; // vector of the keys for the integers representing the edge pairs
    vector<list<morph::Hex>> regionBound; //for each region list of hexes on the boundary
    std::vector<std::vector<std::list<morph::Hex>::iterator>> regionpHexList; //vector of vectors of pHexes (iterators)
    std::map<int,vector<int>> edges; //map of (i,j) edges, uses pair<->int converters
    std::map<int,vector<FLT>> edgeVals; //map of (i,j) edgeVals, uses pair<->int converter
    std::map<int,vector<FLT>> edgeNN; //map of (i,j) edgeNN values, uses pair<->int converter
    vector<vector<FLT> > regionDist; //from each hexdistances to each seed point
    vector<vector<hexGeometry::point>> vCoords; //for each region coordinates of the vertices
    vector<vector<hexGeometry::point>> mCoords; //for each region coordinates of the midpoints
    vector<int> Creg; //for each hex count of different regions it touches
    vector<int> Cnbr; //for each hex count of neighbour hexes
    vector<vector<FLT>> NN, CC; //hold the field values for each hex in a region
    vector<FLT> nn, cc; //hold the field values for each hex in the whole hexGrid
    vector <pair <FLT, FLT>> centres; //seed points for regions
    vector<pair <FLT, FLT>> centroids; // centroids for regions
    vector<std::pair<FLT,FLT>> diff; //difference between seed point and CoG of region
    morph::HexGrid* Hgrid; //original hexGrid
    vector<morph::BezCurvePath<FLT>> curvedBoundary; //vector of boundaries for creating morphed regions
    hexGeometry* hGeo; //supporting class for geometric analysis
    vector<vector<hexGeometry::lineSegment>> radialSegments; //radial segements from original vertices
    vector<vector<FLT>> radialAngles; //angles of the vertices from the centroid.

 /*! class constructor
  * scale controls the resolution of the hexagonal grid
  * xspan controls the width of the hexagonal grid
  * basepath points to the directory for input/output files
  * npoints is the number of Voronoi seed points
  */
    DRegion (int scale, FLT xspan, string basepath, int npoints) {
        this->NUMPOINTS = npoints;
        this->scale = scale;
        this->logpath = basepath;
        this->xspan = xspan;
        FLT s = pow(2.0, scale-1);
        this->ds = 1.0/s;
        this->hexArea = 3.0*tan30*this->ds*this->ds;
        ofstream afile (this->logpath + "/debug.out" );
        ofstream jfile (this->logpath + "/correlateVector.out");
        ifstream bfile(this->logpath + "/centres.inp");
        afile << this->logpath << endl;
        cout << "before creating BezCurve" <<endl;
        srand(time(NULL)); //reseed random number generator
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        this->n = Hgrid->num();
        afile << "after creating HexGrid with " << this->n << " hexes " << endl;
        FLT maxX = Hgrid->getXmax(0.0);
        afile << " the maximum value of x is is " << maxX << endl;
        hGeo = new hexGeometry();
        // now read in the boundary either as a header or as a morph  d
        //   #include "bez5side.h"
        // morph::ReadCurves r("./rat.svg");
        // Hgrid->setBoundary (r.getCorticalPath(),true);
        // this was the original call, I am trying out setBoundaryDregion for debugging
        //    Hgrid->setBoundary(bound,false);
        Hgrid->setEllipticalBoundary (1.0, 1.0);
        cout << "after setting boundary on  H " << Hgrid->num() << endl;
        n = Hgrid->num();
        cout << "after  boundary set HexGrid has " <<  n << " hexes" << endl;
        // now set the centres either read in or randomly generated
        //    #include "centres.h"
        centres.resize(NUMPOINTS);
        centroids.resize(NUMPOINTS);
        for (unsigned int i=0; i<NUMPOINTS; i++) {
            bfile >> centres[i].first >> centres[i].second;
           cout << "centres.first " << centres[i].first << " centres.second " <<  centres[i].second << endl;
        }
        cout << "after setting centres " << endl;
   //these are the vectors of vectors for the regions
        regionDist.resize(n);
        region.resize(n);
        N.resize(n);
        curvedBoundary.resize(NUMPOINTS);
        regionpHexList.resize(NUMPOINTS);
        cout << "after  filleting fish " << " n = " << n <<endl;
        this->Cnbr.resize(n,6); //count of neighbouring hexes
        this->Creg.resize(n,0); //count of neighbouring hexes
        cout << "before hexRegionList"<<endl;
        hexRegionList.resize(n); //neighbouring regions for a hex
        regionList.resize(NUMPOINTS); //neighbouring regions for a region
        regionVertex.resize(NUMPOINTS); //vertices for a region
        radialSegments.resize(NUMPOINTS); // radial line segments
        radialAngles.resize(NUMPOINTS); // radial line segments
        regionBound.resize(NUMPOINTS); //hexes on region boundary unsorted
        sortedBoundary.resize(NUMPOINTS);//hexes on region boundary sorted by angle
        sortedBoundaryNN.resize(NUMPOINTS);//hexes on region boundary sorted by angle
        sortedBoundaryPhi.resize(NUMPOINTS);//hexes on region boundary sorted by angle
        vCoords.resize(NUMPOINTS); //vertex coordinates of a region
        mCoords.resize(NUMPOINTS); //midpoint coordinates of a region
        innerRegion.resize(NUMPOINTS);
      // check the order numbering in hexen
        int hexCount = 0;
        for (auto h : Hgrid->hexen) {
            afile << "loop interation " << hexCount << " h.vi " << h.vi << endl;
            hexCount++;
        }
        // initialise all regions to be innerRegions, this will be corrected later in setCreg()
        for (unsigned int i=0; i<NUMPOINTS; i++) {
            innerRegion[i] = true;
        }

   // making a neighbour array for convenience
        for (int idx = 0; idx < n; idx++) {
            N[idx].resize(6);
            N[idx][0] = Hgrid->d_ne[idx];
            N[idx][1] = Hgrid->d_nne[idx];
            N[idx][2] = Hgrid->d_nnw[idx];
            N[idx][3] = Hgrid->d_nw[idx];
            N[idx][4] = Hgrid->d_nsw[idx];
            N[idx][5] = Hgrid->d_nse[idx];
        }
        cout << "after neighbour array" << endl;

    //arrays to hold the field values
        n = this->Hgrid->num();
        NN.resize(NUMPOINTS);
        CC.resize(NUMPOINTS);
        nn.resize(n);
        cc.resize(n);
        cout << "after alloc nn and CC" <<endl;

        vector <vector <FLT> > sortedDist;
        sortedDist.resize(n);
    // get a list of distances from each Dirichlet point for each hex.
        for (auto h : Hgrid->hexen) {
            for (unsigned int j=0;j<NUMPOINTS;j++) {
                FLT temp = h.distanceFrom(centres[j]);
                regionDist[h.vi].push_back(temp);
            }
        }
       /*!
        *to produce a list of the sorted distances of each hex from the seed points
        *Also produce a vector containing their indices in sorted order
        *Note the process of sorting the indexes in distance order is completely independent
        *from the production of Sorted list but we assume that the similar sorting used for
        *each means that they are compatible.
        */
        for (int i=0;i<n;i++) {
            vector <FLT> tempvector1;
            tempvector1 = regionDist[i];
            std::stable_sort(tempvector1.begin(),tempvector1.end());
            sortedDist[i] = tempvector1;
            vector <int> tempint = sort_indexes(regionDist[i]);
            region[i] = tempint;
        }
        cout << "after region[i] allocated " <<endl;

        /*!
         * Write out the list of hexes that are equidistant from
         * two different nearest seed point.
         */
        for (int i=0;i<n;i++){
            if (sortedDist[i][0] == sortedDist[i][1]){
                afile<<"i ="<<i;
                afile << "  "<<region[i][0]<<" is equidistant from "<< region[i][1];
            }
        }
        cout << "before hex list loop" << endl;

        /*!
         *Determine the outside hexes and set their neighbours to be the central hex
         *if there is no neighbour this is for the boundary conditions of the solvers
         *
         */
        for (auto &h : Hgrid->hexen){
            this->Cnbr.at(h.vi) = 6;
            if (!HAS_nE(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                h.setBoundaryHex();
                N[h.vi][0] = h.vi;
                this->Cnbr.at(h.vi) = -1;
                cout << "no ne" << endl;
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][0]][0]);//nbr region
            }
            if (!HAS_nNE(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                h.setBoundaryHex();
                this->Cnbr.at(h.vi) = -1;
                N[h.vi][1] =  h.vi;
                cout << "no nne" << endl;
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][1]][0]); //nbr region
            }
            if (!HAS_nNW(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                cout << "no nnw" << endl;
                this->Cnbr.at(h.vi) = -1;
                N[h.vi][2] = h.vi;
                h.setBoundaryHex();
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][2]][0]); //nbr region
            }
            if (!HAS_nW(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                h.setBoundaryHex();
                this->Cnbr.at(h.vi) = -1;
                N[h.vi][3] = h.vi;
                cout << "no nw" << endl;
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][3]][0]); //nbr region
            }
            if (!HAS_nSW(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                h.setBoundaryHex();
                this->Cnbr.at(h.vi) = -1;
                N[h.vi][4] = h.vi;
                cout << "no nsw" << endl;
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][4]][0]); //nbr region
            }
            if (!HAS_nSE(h.vi)) {
                hexRegionList[h.vi].push_back(-1);
                h.setBoundaryHex();
                this->Cnbr.at(h.vi)  = -1;
                N[h.vi][5] = h.vi;
                cout << "no nse" << endl;
            }
            else {
                hexRegionList[h.vi].push_back(region[N[h.vi][5]][0]); //nbr region
            }
        } //end of logic determining hex types

           /*!
            * fills the region with appropriate hex indices
            */
       this->regionHex.resize(NUMPOINTS);
       for (unsigned int j=0;j<NUMPOINTS;j++){
           for (auto h : Hgrid->hexen) {
               if (region[h.vi][0] == (int) j) {
                   cout << "hex " << h.vi << " region " << region[h.vi][0] << " j " << j<< endl;
                   this->regionHex[j].push_back(h);
               }
           }
       }


       diff.resize(NUMPOINTS);
       int totalHex=0;
       //print out the regions and count the h/exes in the grid
       for (unsigned int j=0;j<NUMPOINTS;j++){
           totalHex += printRegion(j);
           afile << " total hexes" << totalHex << endl;
       }
       // fill the centroids
       for (unsigned int j=0; j<NUMPOINTS; j++) {
           this->centroids[j] = this->baryCentre(j);
       }
       //set the polar coordinates by region
       for (unsigned int j=0;j<NUMPOINTS;j++){
           diff[j] = this->set_polars(j);
           afile << "diff seed-centre" << diff[j].first << " " << diff[j].second<<endl;
       }
       cout << "after set_polars" << endl;
       //writing hexen as a vector
       for (auto h : Hgrid->hexen) {
           vHexInGrid.push_back(h);
       }
       cout<<"at end of constructor" << " hexes counted " << totalHex << " n " << n << endl;

   } //end of DRegion constructor

   /*!
    * processing of internal boundaries
    */
    void setCreg() {
        for (auto h : this->Hgrid->hexen) {
            int centralRegion = region[h.vi][0];
            int oldRegion = centralRegion;
            int newRegion = 0;
            cout << "just before the Creg logic" << " i " << h.vi << endl;
            for (unsigned int j=0;j<6;j++) {
            //cout << "j = " << j  << " h.vi " << h.vi <<  endl;
                newRegion =  this->hexRegionList[h.vi][j];
                cout << " hexRegionList " << hexRegionList[h.vi][j] << endl;
                if (centralRegion != newRegion) { //its a boundary hex
                    cout << "centralRegion " << centralRegion << " newRegion " << newRegion << endl;
                    h.setBoundaryHex();
                    //this line could have been the problem. Its setting the j neighbour of h.vi to be the central value
                    //I am now setting it to a value that cannot be a h.vi value
                    if (oldRegion != newRegion){ //logic to test if vertex
                        Creg[h.vi]++;
                        cout << " Creg " << Creg[h.vi] << " oldRegion " << oldRegion << " newRegion " << newRegion << endl;
                        oldRegion = newRegion;
                    }
                }
            }
           cout << "end of creg loop" << endl;
        }
    }// end of method setCreg

    int setInnerRegion() {
        int result = 0;
        for (unsigned int i = 0; i<this->NUMPOINTS; i++) {
            unsigned int surroundSize = this->regionList[i].size();
            cout << "surroundSize = " << surroundSize << endl;
            for (unsigned int j=0; j<surroundSize; j++) {
                cout << "in setInnerRegion loop region " << i << " adj region " << regionList[i][j] << endl;
                if (regionList[i][j] == -1) {
                    this->innerRegion[i] = false;
                    result++;
                    break;
                }
            }
        }
        return result;
    }

    void setInternalBoundary() {
        for (auto h : this->Hgrid->hexen) {
            int centralRegion = region[h.vi][0];
            int newRegion = 0;
            cout << "just before the internal boundary logic" << " i " << h.vi << endl;
            for (unsigned int j=0;j<6;j++) {
            //cout << "j = " << j  << " h.vi " << h.vi <<  endl;
                newRegion =  this->hexRegionList[h.vi][j];
                //if (newRegion == -1) continue;
                cout << " hexRegionList " << hexRegionList[h.vi][j] << endl;
                if (centralRegion != newRegion) { //its a boundary hex
                    cout << "centralRegion " << centralRegion << " newRegion " << newRegion << endl;
                    h.setBoundaryHex();
                    //this line could have been the problem. Its setting the j neighbour of h.vi to be the central value
                    //I am now setting it to a value that cannot be a h.vi value
                    N[h.vi][j] = h.vi;
                }
            }
           cout << "end of setInternalBoundary loop" << endl;
        }
    }// end of method setInternalBoundary
    /*!
     * Euclidean distance between two points
     */
    FLT getdist(std::pair<FLT, FLT> a, std::pair<FLT, FLT>  b) {
        FLT result;
        result = sqrt((a.first-b.first)*(a.first-b.first) + (a.second-b.second)*(a.second-b.second));
        return result;
    }

    /*!
    * Finds the indexes from original list in sorted order
    * now using stable_sort to ensure consistency in case of hexes equidistant between
    * centres
    */
    template <typename T>
    vector<int> sort_indexes(const vector<T> &v) {
        vector<int> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        stable_sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2];});
        return idx;
    }

    /*!
     * hashes a pair
     */
    int pair2int (pair<int,int> d, int ibase) {
        int result;
        result = d.first*ibase + d.second;
        return result;
    }

    /*!
     * dehashes the pair
     */
    std::pair<int,int> int2pair(int c, int ibase) {
        std::pair<int, int> result;
        result.first = c / ibase;
        result.second = c%ibase;
        return result;
    }

    /*!
     * test function to ensure hex numbers are not getting scrambled
     */
    void listHex(){
        cout << "in listHex" << endl;
        for (auto h : Hgrid->hexen) {
            cout << "h.vi " << h.vi << endl;
        }
    }

    /*!
     * set the radialSegments for all regions
     */
    void setRadialSegments() {
        hexGeometry::point a,b;
        ofstream zfile (logpath + "/vertexAngles.txt",ios::app);
        for (unsigned int j=0;j<NUMPOINTS;j++) {
            this->radialAngles[j].clear();
            a.first = this->centroids[j].first;
            a.second = this->centroids[j].second;
            unsigned int size = regionVertex[j].size();
            for (unsigned int i=0; i<size; i++) {

                hexGeometry::lineSegment temp;
                b.first = this->Hgrid->d_x[regionVertex[j][i]];
                b.second = this->Hgrid->d_y[regionVertex[j][i]];
                temp = hGeo->createLineSegment(a,b);
                this->vCoords[j].push_back(b);
                radialSegments[j].push_back(temp);
                if (b.second >= a.second)
                {
                     radialAngles[j].push_back(atan2((b.second - a.second) , (b.first - a.first)));
                }
                else
                {
                     radialAngles[j].push_back(2*PI + atan2((b.second - a.second) , (b.first - a.first)));
                }
            }
      // code for all iterations to round corners
           for (unsigned int i=0; i<size; i++)
           {
               FLT x = (vCoords[j][i].first + vCoords[j][(i+1)%size].first)/2.0;
               FLT y = (vCoords[j][i].second + vCoords[j][(i+1)%size].second)/2.0;
               hexGeometry::point p;
               p.first = x;
               p.second = y;
               mCoords[j].push_back(p);
           }
        } //end of loop over regions
        //write out the angles
        for (unsigned int j=0; j<NUMPOINTS; j++) {
            for (unsigned int i=0; i<regionVertex[j].size(); i++) {
                zfile << "vertex angle region " << j << " vertex " << i << " = " << radialAngles[j][i] << endl;
            }
        }
        zfile.close();

    }//end of method setRadialSegments


    /*!
     * Laplacian for solver
     */
    vector<FLT> getLaplacian(vector<FLT> Q, FLT dx) {
        FLT overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<FLT> L(n,0.);
        for(auto h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overds;
        }
        return L;
    }

    /*!
     * Chemotaxis term for Keller-Segel equations
     */
    vector<FLT> chemoTaxis(vector<FLT> Q, vector<FLT> P, FLT dx) {
		vector<FLT> cT(n,0.);
        FLT overds = 1./(1.5*29.0*29.0*dx*dx);

        for (auto h : Hgrid->hexen) {
		    unsigned int i = h.vi;
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
	        FLT dr0Q = (Q[N[i][0]]+Q[i])/2.;
	        FLT dg0Q = (Q[N[i][1]]+Q[i])/2.;
	        FLT db0Q = (Q[N[i][2]]+Q[i])/2.;
	        FLT dr1Q = (Q[N[i][3]]+Q[i])/2.;
	        FLT dg1Q = (Q[N[i][4]]+Q[i])/2.;
	        FLT db1Q = (Q[N[i][5]]+Q[i])/2.;
	  //FLT ncentre = Q[i];

            FLT dr0P = P[N[i][0]]-P[i];
            FLT dg0P = P[N[i][1]]-P[i];
            FLT db0P = P[N[i][2]]-P[i];
	        FLT dr1P = P[N[i][3]]-P[i];
            FLT dg1P = P[N[i][4]]-P[i];
            FLT db1P = P[N[i][5]]-P[i];


	        //finite volume for NdelC, h = s/2
	        cT[i] = (dr0Q*dr0P+dg0Q*dg0P+db0Q*db0P+dr1Q*dr1P+dg1Q*dg1P+db1Q*db1P)*overds;


        } //matches for on i
		return cT;
	} //end of function chemoTaxis

    /*!
     * function to compute the derivative
     */
     void compute_dNNdt(vector<FLT>& inN, vector<FLT>& dNdt, FLT Dn, FLT Dchi) {
         vector<FLT> lapN(this->n,0);
         vector<FLT> cTaxis(this->n,0);
		 //cout << "in compute_dNN just before laplacian" << endl;
        lapN = getLaplacian(inN,this->ds);
        cTaxis = chemoTaxis(inN,cc,this->ds);
        FLT a = 1., b = 1.;
		//cout << "in compute_dNN just before the loop" << endl;
        for (int h=0; h < this->n; h++) {
            dNdt[h] = a-b*inN[h] + Dn*lapN[h] - Dchi*cTaxis[h];
        }
    }

	void compute_dCCdt(vector<FLT>& inC, vector<FLT>&  dCdt, FLT Dc) {
        FLT beta = 5.;
        FLT mu = 1;
		//cout << " before calls to Laplacian " << endl;
        vector<FLT> lapC(this->n,0);
		lapC = getLaplacian(inC,this->ds);
        FLT N2;
        for(int h=0; h < this->n;h++){
            N2 = this->nn[h]*this->nn[h];
            dCdt[h] =  beta*N2/(1.+N2) - mu*inC[h] + Dc*lapC[h];
        }
    }



  //function to timestep coupled equations
    void step(FLT dt, FLT Dn, FLT Dchi, FLT Dc) {
      //dt = dt * 2.5 / Dn;

       // Set up boundary conditions with ghost points
      //cout << " in time step before ghost points" << endl;
 /*     for(auto h : Hgrid->hexen){
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
        if(this->Creg[h.vi] > 0){
          for(int j=0;j<6;j++){
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
             if(N[h.vi][j] < 0) {
       	      // nn[N[h.vi][j]] = nn[h.vi];
                 N[h.vi][j] = h.vi;
			//   cout << " nn " << nn[N[h.vi][j]] << " NN central " << nn[h.vi] << endl;
	           //cc[N[h.vi][j]] = cc[h.vi];
	      }
	     }
	   }
      }
  */
        // 2. Do integration of nn
        {
            // Runge-Kutta integration for A. This time, I'm taking
            // ownership of this code and properly understanding it.

            // Ntst: "A at a test point". Ntst is a temporary estimate for A.
            vector<FLT> Ntst(this->n, 0.0);
            vector<FLT> dNdt(this->n, 0.0);
            vector<FLT> K1(this->n, 0.0);
            vector<FLT> K2(this->n, 0.0);
            vector<FLT> K3(this->n, 0.0);
            vector<FLT> K4(this->n, 0.0);

            /*
             * Stage 1
             */
			 //cout << "in ksSolver before compute_dNNdt" << endl;
            this->compute_dNNdt (this->nn, dNdt, Dn, Dchi);
			 //cout << "in ksSolver after compute_dNNdt" << endl;
            for (int h=0; h< this->n; ++h) {
                K1[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h< this->n; ++h) {
                K2[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dNdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0;h<this->n;h++) {
                this->nn[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
				//this->nn[i] = i * 1.0;
            }
        }

        // 3. Do integration of cc
        {
            // Ctst: "B at a test point". Ctst is a temporary estimate for B.
            vector<FLT> Ctst(this->n, 0.0);
            vector<FLT> dCdt(this->n, 0.0);
            vector<FLT> K1(this->n, 0.0);
            vector<FLT> K2(this->n, 0.0);
            vector<FLT> K3(this->n, 0.0);
            vector<FLT> K4(this->n, 0.0);

            /*
             * Stage 1
             */
            this->compute_dCCdt (this->cc, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K1[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
		    this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K2[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dCdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0; h < this->n; ++h) {
                this->cc[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
            }
        }
        //cout  << "value of nn[5] end Runge " << this->nn[5] <<  " number of hexes " << this->n << endl;

    }//end step


  //function to timestep coupled equations
    void stepOuter(FLT dt, FLT Dn, FLT Dchi, FLT Dc) {
      //dt = dt * 2.5 / Dn;
/*
       // Set up boundary conditions with ghost points
      //cout << " in time step before ghost points" << endl;
      int hcount = 0;
      for(auto &h : Hgrid->hexen){
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
        if(this->Cnbr[h.vi] == -1){
          hcount ++;
          for(int j=0;j<6;j++){
		 //    cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << " Cnbr " << Cnbr[h.vi] << " hcount " << hcount << endl;
             if(N[h.vi][j] == -1) {
       	       nn[N[h.vi][j]] = nn[h.vi];
			//   cout << " nn " << nn[N[h.vi][j]] << " NN central " << nn[h.vi] << endl;
	           cc[N[h.vi][j]] = cc[h.vi];
	      }
	     }
	   }
      }
 */
        // 2. Do integration of nn
        {
            // Runge-Kutta integration for A. This time, I'm taking
            // ownership of this code and properly understanding it.

            // Ntst: "A at a test point". Ntst is a temporary estimate for A.
            vector<FLT> Ntst(this->n, 0.0);
            vector<FLT> dNdt(this->n, 0.0);
            vector<FLT> K1(this->n, 0.0);
            vector<FLT> K2(this->n, 0.0);
            vector<FLT> K3(this->n, 0.0);
            vector<FLT> K4(this->n, 0.0);

            /*
             * Stage 1
             */
			 //cout << "in ksSolver before compute_dNNdt" << endl;
            this->compute_dNNdt (this->nn, dNdt, Dn, Dchi);
			 //cout << "in ksSolver after compute_dNNdt" << endl;
            for (int h=0; h< this->n; ++h) {
                K1[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h< this->n; ++h) {
                K2[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dNdt[h] * dt;
                Ntst[h] = this->nn[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dNdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0;h<this->n;h++) {
                this->nn[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
				//this->nn[i] = i * 1.0;
            }
        }

        // 3. Do integration of cc
        {
            // Ctst: "B at a test point". Ctst is a temporary estimate for B.
            vector<FLT> Ctst(this->n, 0.0);
            vector<FLT> dCdt(this->n, 0.0);
            vector<FLT> K1(this->n, 0.0);
            vector<FLT> K2(this->n, 0.0);
            vector<FLT> K3(this->n, 0.0);
            vector<FLT> K4(this->n, 0.0);

            /*
             * Stage 1
             */
            this->compute_dCCdt (this->cc, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K1[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
		    this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K2[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dCdt[h] * dt;
                Ctst[h] = this->cc[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dCdt[h] * dt;
            }

            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0; h < this->n; ++h) {
                this->cc[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
            }
        }
        //cout  << "value of nn[5] end Runge " << this->nn[5] <<  " number of hexes " << this->n << endl;

    }//end step
    /*!
     * swap the radialSegments for all regions
     * if lvCoords = true take vCoords as the vertices
     * else take mCoords
     */
    void swapRadialSegments(bool lvCoords) {
        hexGeometry::point a,b;
        ofstream zfile (this->logpath + "/vertexAngles.txt",ios::app);
        for (unsigned int j=0;j<NUMPOINTS;j++) {
            this->radialAngles[j].clear();
            a.first = this->centroids[j].first;
            a.second = this->centroids[j].second;
            unsigned int size = regionVertex[j].size();
            for (unsigned int i=0; i<size; i++) {

                hexGeometry::lineSegment temp;
                if (lvCoords) {
                    b = vCoords[j][i];
                }
                else {
                    b = mCoords[j][i];
                }
                temp = hGeo->createLineSegment(a,b);
                radialSegments[j].push_back(temp);
                if (b.second >= a.second)
                {
                     radialAngles[j].push_back(atan2((b.second - a.second) , (b.first - a.first)));
                }
                else
                {
                     radialAngles[j].push_back(2*PI + atan2((b.second - a.second) , (b.first - a.first)));
                }
            }
        } //end of loop over regions
        //write out the angles
        for (unsigned int j=0; j<NUMPOINTS; j++) {
            for (unsigned int i=0; i<regionVertex[j].size(); i++) {
                zfile << "vertex angle region " << j << " vertex " << i << " = " << radialAngles[j][i] << endl;
            }
        }
        zfile.close();
    }

  //function to return area of a region
  FLT regArea (int regNum) {
    FLT area = 0;
    for (unsigned int i=0;i <  this->regionHex[regNum].size();i++){
      area += 1.;
    }
    return area;
  } //end of funtion regArea

  //function to return mean value of NN in the region
  FLT meanNN (int regNum) {
    FLT mean = 0;
    int size = (int) this->regionHex[regNum].size();
    for (int i=0; i<size; i++){
        mean += this->NN[regNum][i];
    }
    if (size != 0) {
        return mean / (1.0 * size);
    }
    else {
         return -999.999;
    }
  } //end of function meanNN

  //function to return mean value of nn in the region
  FLT meannn (int regNum) {
    FLT mean = 0;
    int size = (int) this->regionHex[regNum].size();
    for (auto h : regionHex[regNum]){
        mean += nn[h.vi];
    }
    if (size != 0) {
        return mean / (1.0 * size);
    }
    else {
         return -999.999;
    }
  } //end of function meanNN

  //function to return perimeter of a region
  FLT regPerimeter (int regNum) {
    // cout << "in regPerimeter " << endl;
    FLT perimeter = 0;
    for (auto h : this->regionHex[regNum]){
        if (Creg[h.vi] > 0)
            perimeter += 1.0;
    }
    return perimeter;
  } //end of function regPerimeter


//function to return perimeter of a morphed region
FLT renewRegPerimeter (int regNum) {
  // cout << "in regPerimeter " << endl;
  FLT perimeter = 0;
  for (auto h : this->regionHex[regNum])
    if (h.boundaryHex())
      perimeter += 1.0;
  return perimeter;
} //end of function regPephieter

// method to produes intermediate points in a region boundary
    vector <hexGeometry::point> divideRegionBoundary(int regNum, int ticks) {
        vector <hexGeometry::point> vertices;
        unsigned int size = vCoords[regNum].size();
        cout << "In divideRegionBoundary size " << size << endl;
        for (unsigned int i=0; i < size; i++) {
            FLT xstart = this->vCoords[regNum][i].first;
            FLT ystart = this->vCoords[regNum][i].second;
            FLT xend = this->vCoords[regNum][(i+1)%size].first;
            FLT yend = this->vCoords[regNum][(i+1)%size].second;
            FLT incrX = (xend - xstart)/(1.0*ticks);
            FLT incrY = (yend - ystart)/(1.0*ticks);
            for (int j=0; j<ticks; j++) {
               FLT xval = xstart + j * incrX;
               FLT yval = ystart + j * incrY;
               hexGeometry::point p;
               p.first = xval; p.second = yval;
               vertices.push_back(p);
            }
        }
        cout << " end of  divideRegionBoundary " << endl;
        return vertices;
    }

// method to determine if a point is in a rectangle
    bool inRegion(int regNum, std::pair<FLT, FLT> inPoint, vector<hexGeometry::point> cutter, FLT tol) {
        cout << " in inRegion  region " <<  regNum  << endl;
        bool result;
        unsigned int size = cutter.size();
        cout << "size " << size << " point.x " << inPoint.first << " point.y " << inPoint.second << endl;
        hexGeometry::point testPoint;
        testPoint.first = inPoint.first;
        testPoint.second = inPoint.second;
        vector<FLT> angles;
        FLT windingNumber = 0;
        FLT minAngle = 2 * PI;
        int indexCount = 0;
        for (unsigned int i=0; i<size; i++) {
            hexGeometry::lineSegment firstSide = this->hGeo->createLineSegment(testPoint, cutter[i]);
            hexGeometry::dLine dFirstSide = this->hGeo->segment2dLine(firstSide);
            if (dFirstSide.angle < minAngle) {
                minAngle = dFirstSide.angle;
                indexCount = i;
            }
        }
        cout << " minimum angle " << minAngle << " index " << indexCount;
        cout << " before dLine loop cutter size " << cutter.size() << endl;
        for (unsigned int i=indexCount; i < size + indexCount; i++) {
            hexGeometry::lineSegment tempSide = hGeo->createLineSegment(testPoint, cutter[i%size]);
            hexGeometry::dLine dSide = hGeo->segment2dLine(tempSide);
            FLT correctedAngle = dSide.angle - minAngle;
            angles.push_back(correctedAngle);
            cout << " region " << regNum << " correctedAngle " << correctedAngle << " originalAngle " << dSide.angle << endl;
        }
        cout << " before windingNumber loop" << endl;
        FLT angleSum = 0.0;
        for (unsigned int i=0; i<(size - 1);i++){
            unsigned int lead = (i+1) % size;
            cout << " lead " << angles[lead] << " follow " << angles[i] << " region " << regNum << endl;
            angleSum += angles[lead] - angles[i];
        }
        windingNumber = angleSum / (2.0 * PI);
        if (((1.0 - tol) < windingNumber) && (windingNumber < (1.0 + tol))) {
            result = true;
        }
        else {
            result = false;
        }
        cout << " Winding Number for region " << regNum << " is " << windingNumber << " in region bool " << result << endl;
        return result;
    }


// method to determine if a point is in a polygon
    bool testRegionVertices(int regNum) {
        bool result;
        ofstream outfile ("./logs/testVertices.txt",ios::app);
        unsigned int size = vCoords[regNum].size();
        hexGeometry::point testPoint;
        std::pair<FLT, FLT> inPoint = this->baryCentre(regNum);
        testPoint.first = inPoint.first;
        testPoint.second = inPoint.second;
        outfile << endl;
        outfile << " barycentre " << inPoint.first << " , " << inPoint.second << " vCoords size " << size << endl;
        FLT minAngle = 2 * PI;
        vector<FLT> angles;
        int indexCount = 0;
        for (unsigned int i=0; i<size; i++) {
            hexGeometry::lineSegment firstSide = this->hGeo->createLineSegment(testPoint, vCoords[regNum][i]);
            hexGeometry::dLine dFirstSide = this->hGeo->segment2dLine(firstSide);
            if (dFirstSide.angle < minAngle) {
                minAngle = dFirstSide.angle;
                indexCount = i;
            }
        }
        outfile << " minimum angle " << minAngle << " index " << indexCount;
        outfile << " before dLine loop vertices size " <<  size << endl;
        for (unsigned int i=indexCount; i < size + indexCount; i++) {
            hexGeometry::lineSegment tempSide = hGeo->createLineSegment(testPoint, vCoords[regNum][i%size]);
            hexGeometry::dLine dSide = hGeo->segment2dLine(tempSide);
            FLT correctedAngle = dSide.angle - minAngle;
            angles.push_back(correctedAngle);
            outfile << " i " << i << " correctedAngle " << correctedAngle << " originalAngle " << dSide.angle << endl;
        }
        outfile << " before windingNumber loop" << endl;
        result = true;
        for (unsigned int i=0; i<size;i++){
            unsigned int lead = (i+1) % size;
            if (angles[lead] == 0.0) {
                angles[lead] += 2.0 * PI;
            }
            outfile << " lead " << angles[lead] << " follow " << angles[i] << " region " << regNum << endl;
            if (angles[lead] < angles[i]) {
                result = false;
                break;
            }
        }
        outfile << " region " << regNum << " after windingNumber loopBool  " << result << endl;
        return result;
    }





    // transform vector so its mean is zero
    vector<FLT> meanzero_vector(vector<FLT> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
        vector <FLT> result;
        unsigned int size = invector.size();
        //meanzero << "size " << size << endl;
        FLT sum = 0;
        FLT absSum = 0;
        for (unsigned int i=0; i <size; i++) {
            sum += invector[i];
            absSum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        absSum = absSum / (1.0 * size);
        //meanzero << " mean  " << sum << endl;
        //meanzero << " absolute mean  " << absSum << endl;
        for (unsigned int i=0; i < size; i++) {
            result.push_back(invector[i] - sum);
            //meanzero << " i " << result[i] << endl;
        }
        return result;
    }

    // transform vector by subtracting centval
    vector<FLT> meanzero_vector(vector<FLT> invector, FLT centval) {
        //ofstream meanzero ("meanzero.txt",ios::app);
        vector <FLT> result;
        int size = invector.size();
        //subtract centval from the vector
        for (int i=0; i < size; i++) {
            result.push_back(invector[i] - centval);
            //meanzero << " i " << result[i] << endl;
        }
        return result;
    }

    // transform vector of pairs so its mean is (0,0)
    vector<std::pair<FLT,FLT>> meanzero_vector(vector<std::pair<FLT,FLT>> invector) {
        FLT sum1, sum2 = 0;
        int size = invector.size();
        vector<std::pair<FLT,FLT>> result;
        for (int i=0;i<size;i++) {
            sum1 += invector[i].first;
            sum2 += invector[i].second;
        }
        sum1 = sum1 / (1.0*size);
        sum2 = sum2 / (1.0*size);
        for (int i=0;i<size;i++) {
            std::pair<FLT,FLT> tempPair((invector[i].first - sum1),(invector[i].second - sum2));
            result.push_back(tempPair);
        }
        return result;
    }


    // return the mean of the absolute values of a  vector
        FLT absmean_vector(vector<FLT> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
          FLT sum = 0;
          int size = invector.size();
        //meanzero << "size " << size << endl;
          for (int i=0; i <size; i++) {
            sum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        return sum;
    }

   // return the mean of the values of a  vector
        FLT mean_vector(vector<FLT> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
          FLT sum = 0;
          int size = invector.size();
        //meanzero << "size " << size << endl;
          for (int i=0; i <size; i++) {
            sum += invector[i];
        }
        sum = sum/(1.0*size);
        return sum;
    }

    FLT maxVal( vector<FLT> invector) {
            FLT result = -1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] > result)
                            result = invector[i];
            }
            return result;
    }


    // to find the minimum value of a vector
    FLT minVal( vector<FLT> invector) {
            FLT result = 1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] < result)
                            result = invector[i];
            }
            return result;
    }

//method to normalise a vector
    vector<FLT> normaliseVect( vector<FLT> invector) {
        vector<FLT> result;
        result.resize(0);
        //catch zero size vectors
        if  (invector.size() == 0) {
            cout << "Error in normalise zero size vector" << endl;
            return result;
        }
        FLT vectMean = 0;
        int size = invector.size();
        for(int i=0; i<size; i++) {
            vectMean += invector[i]*invector[i];
        }
        vectMean = sqrt(vectMean / (1.0 * size));
        for (int i=0;i<size;i++) {
            result.push_back(invector[i]/vectMean);
        }
        return result;
    }

//function to print a vector
    void printFLTVect (std::string path, vector<FLT> invect) {
        static int count;
        count ++;
        ofstream Vout (path ,ios::app);
        vector<FLT>::iterator ptr;
        for (ptr = invect.begin(); ptr < invect.end(); ptr++) {
           Vout << *ptr << "  ";
        }
        //Vout << "count " << count << endl;
        //Vout << "---------------------------------------" << endl;
        Vout << endl;
    }


//function to print a vector
    void printIntVect (string path, vector<int> invect) {
        std::ofstream Vout (path ,ios::app);
        vector<int>::iterator ptr;
        for (ptr = invect.begin(); ptr < invect.end(); ptr++) {
           Vout << *ptr << "  ";
        }
        Vout << endl;
    }

  //function to return fraction of area with NN positive
  FLT regNNfrac (int regNum) {
    FLT area = 0;
    FLT positive_area = 0;
    //int size = this->regionHex[regNum].size();
    vector<FLT> normalNN;
    // to normalise the NN field
     for (auto h : this->regionHex[regNum]){
          normalNN.push_back(this->NN[regNum][h.vi]);
      }
      normalNN = this->meanzero_vector(normalNN);
      for (int i=0;i < (int) this->regionHex[regNum].size();i++){
      area += 1.;
      if ((normalNN[i]) > 0) {
          positive_area += 1.0;
      }
    }
    return positive_area / area;
  } //end of function regNNfrac

    int printRegion(int regNum) {
      int result = 0;
      cout << " in printRegion " << regNum << endl;
        for (auto h : this->regionHex[regNum]) {
            cout << h.vi<< " ";
            result++;
            }
        return result;
    }
//function to return fraction of area with nn positive
FLT regnnfrac (int regNum) {
  FLT area = 0;
  FLT positive_area = 0;
  //int size = this->regionHex[regNum].size();
  vector<FLT> normalnn;
  // to normalise the NN field
   for (auto h : this->regionHex[regNum]){
        normalnn.push_back(this->nn[h.vi]);
    }
    normalnn = this->meanzero_vector(normalnn);
    for (int i=0;i < (int) this->regionHex[regNum].size();i++){
    area += 1.;
    if ((normalnn[i]) > 0) {
        positive_area += 1.0;
    }
  }
  return positive_area / area;
} //end of function regnnfrac

// function to give r and theta relative to region centre
    pair <FLT,FLT> set_polars(int regNum){
        pair <FLT, FLT> result;
        result.first = 0.0;
        result.second = 0.0;
        FLT xcentre = this->centroids[regNum].first;
        FLT ycentre = this->centroids[regNum].second;
        FLT minPhi = 10.0;
        FLT maxPhi = -10.0;
        cout << "in set polars region " << regNum << " size " << this->regionHex[regNum].size() << endl;
//set the phi values for each hex, this time relative to the region centre
        for (auto&  h : this->regionHex[regNum]) {
            FLT angle;
            h.r = sqrt((h.x - xcentre)*(h.x -xcentre)
            + (h.y - ycentre)*(h.y - ycentre));
            if ((h.y -ycentre) >= 0) {
              angle =  + atan2((h.y - ycentre), (h.x - xcentre));
              h.phi = angle;
              //cout<< "region" << regNum << " h.phi "  << h.phi<<  " index " << h.vi << endl;
              }
            else {
              angle =  2*PI + atan2((h.y - ycentre), (h.x - xcentre));
              h.phi = angle;
              //cout<< "region " << regNum << " h.phi " << h.phi<<  " index " << h.vi << endl;
              }
            if (angle < minPhi) {
                minPhi = angle;
            }
            if (angle > maxPhi) {
                maxPhi = angle;
            }
        }
        cout << " set_polars max phi " << maxPhi << " minPhi " << minPhi << endl;
        result.first = this->centres[regNum].first - xcentre; //diff between seed point and barycentre
        result.second = this->centres[regNum].second - ycentre ;
        cout << " result.first " << result.first << " result.second " << result.second <<endl;
        return result;
    } //end of function set_polars

    void shift_polars (int regNum, FLT angle) {
        for (auto& h : this->regionHex[regNum]) {
            if (h.phi + angle > 2.0 * PI) {
                h.phi = h.phi + angle - 2*PI;
            }
            else if (h.phi + angle < 0.0) {
                h.phi = h.phi + 2.0 * PI + angle;
            }
            else {
                h.phi = h.phi + angle;
            }
        }
    }


    //function to find all the edges and vertices of the internal boundary
    vector <std::pair<FLT,FLT>> dissectBoundary(void) {
        ofstream hfile ( this->logpath + "/dissectDebug.out" );
        ofstream ifile ( this->logpath + "/regionList.out" );
        ofstream kfile ( this->logpath + "/edgesList.out" );
        ofstream lfile ( this->logpath + "/verticesList.out" );
        ofstream ufile ( this->logpath + "/keysList.out" );
        hfile<<"just in dissectBoundary"<<endl;
        vector<std::pair<FLT,FLT>> result;
        vector<int>  regionBoundary; //holding array for the indices of hexes in each region boundary
        this->edgeIndex.resize(0); //reset edgeIndex
        int sideCount = 0;
        for (unsigned int iregion=0; iregion < NUMPOINTS; iregion++) { //loop over regions
            cout << " region index " << iregion << " size "  <<this->regionHex[iregion].size() << endl;
            vector<FLT> rB; //holds the angles of the hex from the centroid
            rB.resize(0);
            regionBoundary.resize(0);
// fill the regionBoundary and record polar angles in rB
            for (auto h : this->regionHex[iregion]) {
                FLT angle;
                if (this->Creg[h.vi] >0){
                    regionBoundary.push_back(h.vi);
                    angle = h.phi;
                    cout<< " getPhi test " << angle <<  " index " << h.vi << endl;
                    rB.push_back(angle);
                }
            } //end of loop on a single region

            cout<<"after filling of regionBoundary" <<endl;
            cout<<"regionBoundary.size "<<regionBoundary.size() << endl;
            vector<int> irB; //holds the sorted boundary indicies
            irB = sort_indexes(rB); //irB holds sorted indices after sort on theta
            for (unsigned int idx=0; idx<irB.size();idx++) {
                this->sortedBoundary[iregion].push_back(regionBoundary[irB[idx]]);
            }
            hfile << " irB size " << irB.size() << " rB size " << rB.size() << endl;
//            for (unsigned int i=0;i<irB.size();i++){
//                hfile << " Creg " << Creg[regionBoundary[irB[i]]] << " index " << regionBoundary[irB[i]] << " theta " << rB[irB[i]] <<endl;
//            } //debugging loop for dissectdebug
            unsigned int irBSize = irB.size();
            if (irBSize == 0) {
                cout << "WARNING: region i " << region[iregion][0] << " irB size " << irBSize << endl;
                continue;
            } // catches empty boundaries.

// we now walk round the boundary in theta order, find the first vertex then proceed from there to find vertices in order
            unsigned int idissect = 0; //counts number of boundary hexes processed
            unsigned int Vcount = 0; //count of the vertices
            unsigned int Ecount = 0; //count of the edges
            int newVertex; //integer to enumerate the vertices
            unsigned int offset = 0; //number of boundary hexes before the first vertex
            vector<int> ihE; //contains the sorted indicies of each edge,
                //find the offset to the first vertex
            while (offset < irBSize) {
                if (Creg[regionBoundary[irB[offset]]] > 1)
                    {
                        //Vcount++; //its a vertex
                        newVertex = regionBoundary[irB[offset]];
                        idissect++;
                        break; // found the first vertex
                    }
                offset++;
            }
            cout<<"after offset loop" << " offset " << offset  << " idissect " << idissect << endl;
            while ((idissect < irBSize)) {
                Ecount = 0;
                ihE.resize(0);
    // while loop to catch the nasty case of artificial adjacent vertices, we only count the last.
    // this is a result of our proceeding via hex body rather than hex edge.

                while (Creg[regionBoundary[irB[(idissect + offset) % irBSize]]] > 1) {
                    newVertex = regionBoundary[irB[(idissect+offset)%irBSize]];
                    cout << "in vertex loop" << " idissect " << idissect << " Creg "
                    << Creg[regionBoundary[irB[(idissect + offset) % irBSize]]] << endl;
                    idissect++;
                    Vcount++; //count all vertices even if adjacent
                } //end of loop to trap adjacent vertices
    //walk along the edge until the next vertex
                cout << "Creg " << Creg[regionBoundary[irB[(idissect + offset)%irBSize]]] << " boundary " << regionBoundary[irB[(idissect + offset)%irBSize]] << " newVertex Creg " << Creg[newVertex] << endl;
                regionVertex[iregion].push_back(newVertex);
                //now fill the edge until the next vertex is encountered.
                while ((this->Creg[regionBoundary[irB[(idissect + offset)%irBSize]]] == 1) && (idissect < irBSize)) {
                    ihE.push_back(regionBoundary[irB[(idissect + offset)%irBSize]]);
                    Ecount++;
                    idissect++;
                }
                ihE.insert(ihE.begin(),newVertex); //edge contrains the start vertex but not the end vertex
                Ecount++;
                cout << "after edge loop Ecount " << Ecount << " Vcount " << Vcount <<endl;
                if (Ecount == 0) {
                    cout<<"WARNING - empty edge=========================================="<<endl;
                    continue;
                } //loop to catch empty edge
                cout<<"ihE Size "<<ihE.size() <<endl;
                cout <<"Vcount "<< Vcount << " Ecount "<< Ecount << endl;
                int regMiddle = region[ihE.rbegin()[1]][0]; // region of the penultimate hex in the edge
                if (regMiddle != (int) iregion) {
                    cout << "ERROR: penultimate hex in edge has different region" << endl;
                }
                int edgeOuter = -2;
                //find the first region not the same as the central, since Creg = 1 there can only be one such region.
                for (int ihex = 0; ihex<6; ihex++){
                    if (regMiddle != hexRegionList[ihE.rbegin()[1]][ihex]){
                        edgeOuter = hexRegionList[ihE.rbegin()[1]][ihex];
                        break;
                    }
                }
                hfile<<"after edgeOuter assignment"<<endl;
                hfile<<"edgeOuter "<<edgeOuter<< " edgeInner " << iregion << endl;
                if (edgeOuter > -1) { // edgouter = -1 means that the edge is on the outside of the computational region
                    std::pair <int,int> keypair(iregion,edgeOuter);
                    int keyint = this->pair2int(keypair,this->base);
                    std::pair <int, vector<int>> p1(keyint,ihE);
                    hfile << "region " << iregion << " after pair set keyint " << keyint << " edgeOuter " << edgeOuter << endl;
                    if (this->edges.insert(p1).second) {
                        this->edgeIndex.push_back(keyint);
                    }
                    sideCount++;
                    hfile << "after edges insert edge size " << ihE.size() << endl;
                    hfile <<"=================================================="<<endl;
                    this->regionList[iregion].push_back(edgeOuter);
                }
                else {
                    hfile << "edgeouter "<< edgeOuter << " edge size " << ihE.size() << endl;
                    this->regionList[iregion].push_back(edgeOuter);
                }
            } // end of idissect loop
            hfile << " after idissect loop region " << iregion<< " regionList size " << regionList[iregion].size()<<  endl;
            hfile << "*********************************************************************************************" << endl;
            if (regionList[iregion].size() != regionVertex[iregion].size()) {
                hfile << "WARNING: duplicate region in regionList for region " << iregion <<  " Vcount " << Vcount << endl;
                hfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            }
            if (regionList[iregion].size() == 0) {
                hfile << "WARNING: region  " << iregion << " has no neighbouring regions" << endl;
                hfile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                continue;
            }
    //write out to  regionList file the neighbouring regions to this one
            ifile << "number of nbrs for region " << iregion << " is " << regionList[iregion].size() << endl;
            for (unsigned int inbr = 0; inbr < regionList[iregion].size(); inbr++) {
                ifile << " r " << iregion << " rNbr " << regionList[iregion][inbr];
                cout << " r " << iregion << " rNbr " << regionList[iregion][inbr];
                ifile << endl;
                ifile << "---------------------------------------"<< endl;
            }
            cout <<endl;
    // write to vertexlist file, vertex x and y coordinates
            lfile << "number of vertices for region " << iregion << " is " << regionVertex[iregion].size() << endl;
            for (unsigned int idx=0; idx < irB.size(); idx++){
                for (unsigned int ivtx=0; ivtx < regionVertex[iregion].size(); ivtx++) {
                    cout << " in vtx loop " << ivtx << " idx " << idx << " vertex " << regionVertex[iregion][ivtx] << " irB " << irB[idx] << endl;
                    if (regionBoundary[irB[idx]] == regionVertex[iregion][ivtx]) {
                        lfile << " vertex " << regionVertex[iregion][ivtx] << " x " << Hgrid->d_x[regionVertex[iregion][ivtx]] << " y " << Hgrid->d_y[regionVertex[iregion][ivtx]] << " theta " << rB[irB[idx]];
                        lfile << " vertex regions " << region[regionVertex[iregion][ivtx]][0] << " , " <<  region[regionVertex[iregion][ivtx]][1] <<  " , " << region[regionVertex[iregion][ivtx]][2] <<endl;
                    }
                }
            }
            lfile << endl;
            lfile << "---------------------------------------"<< endl;
        } //end of loop on regions
        cout << " after loop on regions edgeIndex size " << edgeIndex.size() << " edges size " << edges.size() << " sideCount " << sideCount << endl;
        int difference = 0;
// loops to print out the edges map structure
        int countIndex = 0;
        int halfway = int(edgeIndex.size()) / 2;
        for (unsigned int irlook = 0; irlook < NUMPOINTS; irlook++) {
            for (unsigned int jrlook = 0; jrlook < NUMPOINTS; jrlook++) {
                std::pair <int,int> klook(irlook,jrlook);
                std::pair <int, int> koolk(jrlook,irlook);
                int k = pair2int(klook,this->base);
                int k1 = pair2int(koolk,this->base);
                if (edges.count(k) == 0) // no entry in the container for K
                    continue;
                else {
                    int sizeij = edges[k].size();
                    int sizeji = edges[k1].size();
                    difference = (sizeij - sizeji)*(sizeij - sizeji);
                    if ((sizeij*sizeji != 0) && (difference < 100000000)) { // just to check that the edge size hasnt run away
                        kfile << irlook << " " << jrlook << " sizeij " << sizeij << " "<< jrlook << " " << irlook << " sizeji " << sizeji << endl;
                        hfile << " k = " << k << " k1 = " << k1 << " countIndex = " << countIndex << " edge Indexij = " << edgeIndex[countIndex] << " edgeIndexji " <<  edgeIndex[countIndex + halfway] <<  endl;
                        countIndex++;
                    }
                    else {
                        hfile << " size was zero for k " << k << " and k1 " << k1 << " countIndex " << countIndex << " edge Indexij = " << edgeIndex[countIndex] << " edgeIndexji " <<  edgeIndex[countIndex + halfway]  << endl;
                        countIndex++;
                    } //end of condition edge non empty
                } //end of condition no edge at this index
            } // end of inner loop over edges
        } //end of outer loop over edges
        for (unsigned int i = 0; i<edgeIndex.size(); i++) {
            ufile << " for number in edgeIndex " << i << " value is " << edgeIndex[i] << endl;
        }
        return result;
    }//end of function dissectBoundary

//method to identify vertices of the Voronoi tessellation
    void voronoiVertices(int region)
    {
        ofstream vtxFile(this->logpath + "/voronoiVertices.out",ios::app);
        vtxFile << "number of vertices for region " << region << " is " << this->regionVertex[region].size() << endl;
        for (unsigned int ivtx=0; ivtx < this->regionVertex[region].size(); ivtx++) {
            vtxFile << " vertex " << this->regionVertex[region][ivtx] << " x " << this->Hgrid->d_x[this->regionVertex[region][ivtx]] << " y " << this->Hgrid->d_y[this->regionVertex[region][ivtx]];
            vtxFile << " vertex regions " << this->region[this->regionVertex[region][ivtx]][0] << " , " <<  this->region[this->regionVertex[region][ivtx]][1] <<  " , " << this->region[this->regionVertex[region][ivtx]][2] <<endl;
        }
        vtxFile << endl;
    } // end of method voronoiVertices




    /*
     * function to correlate matching edges
     */
    FLT correlate_edges(int iter, bool lzeroMean=true)
    {
        FLT result = 0;
        string str = to_string(iter);
        ofstream edgefile(this->logpath + "/edgeCorrelations" + str + ".txt", ios::app);
        ofstream edgerest(this->logpath + "/edgeRest" + str + ".txt");
        string dataname = "/correlate" + str + ".data";
        ofstream correl(this->logpath + dataname,ios::app);

        vector<FLT> first;
        vector<FLT> second;
        vector<FLT> dinterp;
        edgefile << " In correlate_edges " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        for (unsigned int i = 0; i <NUMPOINTS; i++)
        {

            FLT nnmean1, nnmean2;
            nnmean1 = this->meannn(i); //find the mean of nn in the region
            edgefile << " mean of nn in region " << i << "is " << nnmean1 << endl;
            for (auto j = this->regionList[i].begin(); j < this->regionList[i].end();j++)
            {
                if (*j == -1) continue;
                nnmean2 = this->meannn(*j); //find the mean of nn in the region
                edgefile << " j iteration " << *j << " nnmean2 " << nnmean2 << endl;
                first.resize(0);
                second.resize(0);
                dinterp.resize(0);
                std::pair<int,int> edgePair1(i,*j);
                int edgeIndex1 = this->pair2int(edgePair1,this->base);
                std::pair<int,int>  edgePair2(*j,i);
                int edgeIndex2 = this->pair2int(edgePair2,this->base);
                edgefile << "region " << i << " outer " << *j << " index1 " << edgeIndex1 << " index2 " << edgeIndex2 << endl;
                int count1 = 0;
                int count2 = 0;
                FLT ratio;
                FLT correlationValue = 0;
                // first and 2 have the values of nn on the edge with the mean of the region subtracted
                //create two vectors from the opposing edges
                for (auto itr = this->edges[edgeIndex1].begin(); itr != this->edges[edgeIndex1].end();itr++)
                {
                    first.push_back(this->nn[*itr]);
                    count1++;
                }
                first = this->meanzero_vector(first, nnmean1);
                for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++)
                {
                    second.push_back(this->nn[*itr]);
                    count2++;
                }
                second = this->meanzero_vector(second, nnmean2);
                //vectors are indexed in opposite directions either side
                std::reverse(second.begin(),second.end());
                for (unsigned int  i = 0; i < first.size(); i++) {
                    if (isnan(first[i]))
                        edgefile << "Nan at " << i << " in first" << endl;
                }
                for (unsigned int  i = 0; i < second.size(); i++) {
                    if (isnan(second [i]))
                        edgefile << "Nan at " << i << " in second" << endl;
                }
                //Do this if the vectors are the same size
                if (first.size() == second.size() && first.size()*second.size() != 0)
                {
                    correlationValue = this->correlate_Eqvector(first, second, lzeroMean);
                    ratio = 1.0;
                    edgefile << " if 1 region " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
                    //edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << endl;
                } //end of code if both edges are equal
                //do this if the first vector is bigger than the second
                else if (first.size() > second.size() && first.size()*second.size() != 0)
                {
                    dinterp = this->equalize_vector(second,first);
                    correlationValue = this->correlate_Eqvector(first, dinterp, lzeroMean);
                    ratio = 1.0 * second.size() / (1.0 * first.size());
                    if (correlationValue > -2) {
                        result += fabs(correlationValue);
                        countResult++;
                        correl << correlationValue << "  " <<  endl;
                    }
                    edgefile << " if 2 region " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
                }
                //do this if the first vector is smaller than the second.
                else if (first.size() < second.size() && first.size()*second.size() != 0)
                {
                    dinterp = this->equalize_vector(first,second);
                    correlationValue = this->correlate_Eqvector(dinterp, second, lzeroMean);
                    ratio = 1.0 * first.size() / (1.0 * second.size());
                    if (correlationValue > -2) {
                        result += fabs(correlationValue);
                        countResult++;
                        correl << correlationValue << "  " <<  endl;
                    }
                    edgefile << " if 3 region " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
                }
                //something has gone wrong
                else
                {
                    edgefile << " if 4 region " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
                }
            } //end of single edge comparison
        } // end of loop on regions
        edgefile << " countResult "<<countResult<< " entries in edgeIndex " << this->edgeIndex.size() << " entries in edges " << this->edges.size() << endl << endl;;
        result = result / (countResult * 1.0);
        edgefile.close();
        correl.close();
        return result;
    } //end of function correlate_edges
/*
 * function to insert cosine functions along the edges of the tesselation
 * phaseShift controls how the edges are shifted relative to each other
 * mode alters the frequency of the cosine function.
 */
    void insert_cosines(FLT phaseShift, FLT mode, bool Ljiggle = false) {
        string edgeVal = logpath + "/edgeVals.out";
        this->edgeVals.clear();
        unsigned int seed = time(NULL);
        morph::RandUniform<FLT> ruf(seed);
        for (unsigned int i = 0; i <NUMPOINTS; i++) {
            for (auto j = this->regionList[i].begin(); j != this->regionList[i].end();j++) {
                int count1 = 0;
                std::pair<int,int> edgePair(i,*j);
                int edgeIndex = this->pair2int(edgePair,this->base);
                int s1 = this->edges[edgeIndex].size();
// tempvect1 and 2 have the values of NN on the edge with the mean of the region subtracted
                FLT xstep = (PI * 1.0) / (1.0 * (s1 - 1));
                FLT xval = 0;
                int xcount = 0;
                FLT amode = 0;
                vector<FLT> tempvect;
                FLT jump;
                FLT hold = ruf.get();
                //cout << "hold = " << hold << endl;
                if (hold < 0.5)
                    jump = 1.0;
                else
                    jump = 0.0;
                amode = mode * jump + 1.0;
                //amode = mode + 1.0;
                cout << "amode " << amode << endl;
                FLT phase;
                if (ruf.get() <= 0.5)
                    phase =   phaseShift * PI * ruf.get() ;
                else
                    phase = - phaseShift * PI * ruf.get();
                //for (auto itr = this->edges[edgeIndex].begin(); itr < this->edges[edgeIndex].end();itr++) {
                for (int i=0; i< s1; i++) {
                    xval = xcount * xstep + phase;
                    FLT val = cos (amode * xval);
                    tempvect.push_back(val);
                    count1++;
                    xcount++;
                }
                xcount--;
                if (ruf.get() < 0.5) {
                    std::reverse(tempvect.begin(),tempvect.end()); //edge is oriented reverse direction
                }
                std::pair <int, vector<FLT>> p(edgeIndex,tempvect);
                this->edgeVals.insert(p);

                printFLTVect(edgeVal,tempvect);
                cout << "size of edgeVal region " << i << " edge " << *j << " is " << tempvect.size() << " xstep is " << xstep << " pi " << xstep*xcount << endl;
            } //end of loop over edges of a region
        } // end of loop over all regions
        cout << "in insert_cosines edges size " << this->edges.size() << " edgeVals size " << this->edgeVals.size() << endl;
    } // end of function insert cosines

/*
 * function to correlate adjacent edges
 */
    FLT adjacent_cosines() {
        FLT result = 0;
        ofstream edgefile(this->logpath + "/edgeCorrelations0.txt",ios::app);
        ofstream edgerest(this->logpath + "/edgeRest.txt");
        ofstream correl(this->logpath + "/correlate0.data",ios::app);
        vector<FLT> first;
        vector<FLT> second;
        vector<FLT> dinterp;
        edgefile << " In adjacent cosines " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        int printInt = 15;
        std::string filei = logpath + "/ival.Vect";
        std::string filej = logpath + "/jval.Vect";
        ofstream iout (filei,ios::app);
        ofstream jout (filej,ios::app);
        unsigned int seed = time(NULL);
        // A rando2yym uniform generator returning real/FLTing point types
        FLT correlationValue = 0;
        morph::RandUniform<FLT> ruf(seed);
        for (unsigned int i = 0; i <NUMPOINTS; i++) {
            for (auto j = this->regionList[i].begin(); j != this->regionList[i].end();j++) {
                //edgefile << " j iteration " << *j << endl;
                first.resize(0);
                second.resize(0);
                dinterp.resize(0);
                std::pair<int,int> edgePair1(i,*j);
                int edgeIndex1 = pair2int(edgePair1, this->base);
                first = this->edgeVals[edgeIndex1];
                std::pair<int,int>  edgePair2(*j,i);
                int edgeIndex2 = pair2int(edgePair2, this->base);
                second = this->edgeVals[edgeIndex2];
                int s1 = first.size();
                int s2 = second.size();
                int s3 = 0;
                if ((s1 < 5) || (s2 < 5)) {
                    continue;
                }
// tempvect1 and 2 have the values of NN on the edge with the mean of the region subtracted
                if (s1 == s2 && s1*s2 != 0)
                {
                    correlationValue = this->correlate_Eqvector(first, second);
                    result += fabs(correlationValue);
                    countResult++;
                    correl << correlationValue << endl;
                    edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << " cV " << correlationValue <<endl;
                } //end of code if both edges are equal
                else if (s1 > s2 && s1*s2 != 0)
                {
                    dinterp = this->equalize_vector(second, first);
                    s3 = dinterp.size();
                    if (s3 == 0) {
                        edgefile << " i " << i << " count1 " << first.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    }
                    correlationValue = this->correlate_Eqvector(dinterp,first);
                    if (correlationValue > -2) {
                        result += fabs(correlationValue);
                    if (countResult%printInt == 0) {
                        iout << " Correlation Value = " << correlationValue << endl;
                        //printFLTVect(filei, first);
                        jout << " Correlation Value = " << correlationValue << endl;
                        //printFLTVect(filej, dinterp);
                    }
                        countResult++;
                        edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << " cV " << correlationValue <<endl;
                        correl << correlationValue << endl;
                    }
                    else {
                        edgefile << "ERROR: edges " << i << " and " << *j << " have not been equalised" << endl;
                    }

                }
                else if (s1 < s2 && s1*s2 != 0)
                {
                    dinterp = this->equalize_vector(first, second);
                    s3 = dinterp.size();
                    if (s3 == 0){
                        edgefile << " i " << i << " count2 " << second.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    }
                    correlationValue = this->correlate_Eqvector(dinterp, second);
                    if (correlationValue > -2) {
                        result += fabs(correlationValue);
                    if (countResult%printInt == 0) {
                        iout << " Correlation Value = " << correlationValue << endl;
                        //printFLTVect(filei, dinterp);
                        jout << " Correlation Value = " << correlationValue << endl;
                        //printFLTVect(filej, second);
                    }
                        countResult++;
                        edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << " cV " << correlationValue <<endl;
                        correl << correlationValue << endl;
                    }
                    else {
                        edgefile << "ERROR: edges " << i << " and " << *j << " have not been equalised" << endl;
                    }
                }
                else
                {
                    edgefile << "error zero size for one of the edges " << i << " size " << first.size() << " second " << *j << " size " << second.size() << endl;
                }
            } //end of single edge comparison
        } // end of loop on regions
        edgefile << " countResult "<<countResult<<endl;
        result = result / (countResult * 1.0);
        edgefile.close();
        return result;
    } //end of function adjacent_cosines


    // method to compare random pairs of edges
    void random_correlate(const int max_comp, const int morphNum, bool lNN=true, bool lzero=true) {
        ofstream jfile;
        ofstream kfile;
        int max_rand = edgeIndex.size();
        string str = to_string(morphNum);
        jfile.open(this->logpath + "/random_correlate" + str + ".txt",ios::app);
        kfile.open(this->logpath + "/random_correlate" + str + ".data",ios::app);
        vector<FLT> dinterp, first, second;
        FLT corr;
        jfile << " max_rand " << max_rand << " max_comp " << max_comp << endl;
        /*
         * extract the region and neighbour region of the edge
         */
        int count = 0;
        while (count <max_comp) {
           int r1 = rand() % max_rand;
           int r2 = rand() % max_rand;
           int rr1 = this->edgeIndex[r1]; //edgeIndex is a vector of the integer keys of edges
           int rr2 = this->edgeIndex[r2];
           jfile << " r1 " << r1  << " r2 " << r2 << " rr1 " << rr1 << " rr2 " << rr2 << endl;
           int s3 = 0;
           first.resize(0);
           second.resize(0);
           int s1 = edges[rr1].size(); //edge value is integer array of the hex identifiers of the edge
           int s2 = edges[rr2].size();
           int reg1 = rr1 / 1000;
           int reg2 = rr2 / 1000;
           int out1 = rr1 % 1000;
           int out2 = rr2 % 1000;
           jfile << "In random_correlate region 1 " << reg1 << " rr1 " << rr1 << " region 2 " << reg2 << " rr2 " << rr2 << endl;
           //check if edges are from the same region or from regions that are adjacent
           if ((reg1 == reg2) || (reg1 == out2) || (reg2 == out1)) {
              jfile << "in random_correlate neighbour detected reg 1 " << reg1 << " reg 2  " << reg2 << endl;
              continue;
           }
           //need the mean to normalise the boundary vectors
           //FLT NNmean1 = this->meanNN(reg1);
           FLT NNmean1=0;
           cout << "after meanNN call" << endl;
           //FLT NNmean2 = this->meanNN(reg2);
           FLT NNmean2=0;
           if ((s1 != 0) && (s2 != 0)) //neither edge is empty
           {
               if (lNN) {  //we are working on the NN field as a vector<vector>>
                   jfile  << "rr1 " << rr1 << " s1 " << s1 << " rr2 " << rr2 <<  " s2 " << s2 << endl;
                   for (auto itr = this->edges[rr1].begin(); itr != this->edges[rr1].end();itr++) {
                       first.push_back(this->NN[reg1][*itr]);
                   }
                   first = this->meanzero_vector(first, NNmean1);
                   for (auto itr = this->edges[rr2].begin(); itr != this->edges[rr2].end();itr++) {
                       second.push_back(this->NN[reg2][*itr]);
                   }
                   second = this->meanzero_vector(second,NNmean2);
               }

               else{  //we are working on the nn field as a vector
                   jfile  << "rr1 " << rr1 << " s1 " << s1 << " rr2 " << rr2 <<  " s2 " << s2 << endl;
                   for (auto itr = this->edges[rr1].begin(); itr != this->edges[rr1].end();itr++) {
                       first.push_back(this->nn[*itr]);
                   }
                   first = this->meanzero_vector(first, NNmean1);
                   for (auto itr = this->edges[rr2].begin(); itr != this->edges[rr2].end();itr++) {
                       second.push_back(this->nn[*itr]);
                   }
                   second = this->meanzero_vector(second,NNmean2);
               }
               if (s1 < s2) {
                   dinterp = equalize_vector(first , second);
                   s3 = dinterp.size();
                   corr = correlate_Eqvector(dinterp, second, lzero);
                   if (corr == -2) {
                       continue;
                   }
               }
               else if (s1 > s2) {
                   dinterp = equalize_vector(second, first);
                   s3 = dinterp.size();
                   corr = correlate_Eqvector(dinterp, first, lzero);
                   if (corr == -2) {
                       continue;
                   }
               }
               else {
                   corr = correlate_Eqvector(first, second, true);
                   s3 = 0;
               }
               jfile <<  " r1 " << r1 << " rr1 " << rr1 <<" s1 " << s1 << " r2 " << r2 << " rr2 " << rr2 << " s2 " << s2 << " s3 " << s3 << " correlate " << corr << endl << endl;
               kfile << corr << endl;
           } // end of if test for empty edge
        count++;
        } // end of while loop
    }


    // method to compare random pairs of edges for Special functions
    void Srandom_correlate(const int max_comp, const int morphNum, bool lNN=true, bool lzero=true) {
        ofstream jfile;
        ofstream kfile;
        int max_rand = edgeIndex.size();
        string str = to_string(morphNum);
        jfile.open(this->logpath + "/random_correlate" + str + ".txt",ios::app);
        kfile.open(this->logpath + "/random_correlate" + str + ".data",ios::app);
        vector<FLT> dinterp, first, second;
        FLT corr;
        jfile << " max_rand " << max_rand << " max_comp " << max_comp << endl;
        /*
         * extract the region and neighbour region of the edge
         */
        int count = 0;
        while (count <max_comp) {
           int r1 = rand() % max_rand;
           int r2 = rand() % max_rand;
           int rr1 = this->edgeIndex[r1]; //edgeIndex is a vector of the integer keys of edges
           int rr2 = this->edgeIndex[r2];
           jfile << " r1 " << r1  << " r2 " << r2 << " rr1 " << rr1 << " rr2 " << rr2 << endl;
           int s3 = 0;
           first.resize(0);
           second.resize(0);
           int s1 = edges[rr1].size(); //edge value is integer array of the hex identifiers of the edge
           int s2 = edges[rr2].size();
           int reg1 = rr1 / 1000;
           int reg2 = rr2 / 1000;
           int out1 = rr1 % 1000;
           int out2 = rr2 % 1000;
           jfile << "In random_correlate region 1 " << reg1 << " rr1 " << rr1 << " region 2 " << reg2 << " rr2 " << rr2 << endl;
           //check if edges are from the same region or from regions that are adjacent
           if ((reg1 == reg2) || (reg1 == out2) || (reg2 == out1)) {
              jfile << "in random_correlate neighbour detected reg 1 " << reg1 << " reg 2  " << reg2 << endl;
              continue;
           }
           //need the mean to normalise the boundary vectors
           //FLT NNmean1 = this->meanNN(reg1);
           FLT NNmean1=0;
           cout << "after meanNN call" << endl;
           //FLT NNmean2 = this->meanNN(reg2);
           FLT NNmean2=0;
           int edge1Size =edges[rr1].size();
           int edge2Size =edges[rr2].size();
           if ((s1 != 0) && (s2 != 0)) //neither edge is empty
           {
               if (lNN) {  //we are working on the NN field as a vector<vector>>
                   jfile  << "rr1 " << rr1 << " s1 " << s1 << " rr2 " << rr2 <<  " s2 " << s2 << endl;
                   for (int i=0; i<edge1Size; i++) {
                       first.push_back(this->NN[reg1][i]);
                   }


                   first = this->meanzero_vector(first, NNmean1);
                   for (int i=0; i<edge2Size; i++) {
                       second.push_back(this->NN[reg2][i]);
                   }
                   second = this->meanzero_vector(second,NNmean2);
               }

               else{  //we are working on the nn field as a vector
                   jfile  << "rr1 " << rr1 << " s1 " << s1 << " rr2 " << rr2 <<  " s2 " << s2 << endl;
                   for (auto itr = this->edges[rr1].begin(); itr != this->edges[rr1].end();itr++) {
                       first.push_back(this->nn[*itr]);
                   }
                   first = this->meanzero_vector(first, NNmean1);
                   for (auto itr = this->edges[rr2].begin(); itr != this->edges[rr2].end();itr++) {
                       second.push_back(this->nn[*itr]);
                   }
                   second = this->meanzero_vector(second,NNmean2);
               }
               if (s1 < s2) {
                   dinterp = equalize_vector(first , second);
                   s3 = dinterp.size();
                   corr = correlate_Eqvector(dinterp, second, lzero);
                   if (corr == -2) {
                       continue;
                   }
               }
               else if (s1 > s2) {
                   dinterp = equalize_vector(second, first);
                   s3 = dinterp.size();
                   corr = correlate_Eqvector(dinterp, first, lzero);
                   if (corr == -2) {
                       continue;
                   }
               }
               else {
                   corr = correlate_Eqvector(first, second, true);
                   s3 = 0;
               }
               jfile <<  " r1 " << r1 << " rr1 " << rr1 <<" s1 " << s1 << " r2 " << r2 << " rr2 " << rr2 << " s2 " << s2 << " s3 " << s3 << " correlate " << corr << endl << endl;
               kfile << corr << endl;
           } // end of if test for empty edge
        count++;
        } // end of while loop
    }
    // method to compare random pairs of edges
    void random_cosines(const int max_comp, const int morphNum) {
        ofstream jfile;
        ofstream kfile;
        int max_rand = edgeIndex.size();
        max_rand = 205;
        string str = to_string(morphNum);
        jfile.open(this->logpath + "/random_correlate" + str + ".txt",ios::app);
        kfile.open(this->logpath + "/random_correlate" + str + ".data",ios::app);
        string vectfile = logpath + "/random.Vect";
        vector<FLT> dinterp, first, second;
        FLT corr;
        unsigned int seed = time(NULL);
        // A rando2yym uniform generator returning real/FLTing point types
        morph::RandUniform<FLT> ruf(seed);
        jfile << " max_rand " << max_rand << " max_comp " << max_comp << endl;
        /*
         * extract the region and neighbour region of the edge
         */
        int count = 0;
        while (count <max_comp) {
           int r1 = rand() % max_rand;
           int r2 = rand() % max_rand;
           jfile << "jsst before interger assingnments" << endl;
           int rr1 = this->edgeIndex[r1]; //edgeIndex is a vector of the integer keys of edges
           int rr2 = this->edgeIndex[r2];
           jfile << " r1 " << r1  << " r2 " << r2 << " rr1 " << rr1 << " rr2 " << rr2 << endl;
           int s3 = 0;
           first.resize(0);
           second.resize(0);
           jfile << "jsst before size assingnments" << endl;
           int s1 = this->edgeVals[rr1].size(); //edge value is integer array of the hex identifiers of the edge
           int s2 = this->edgeVals[rr2].size();
           if ((s1<5) || (s2<5)) {
               continue;
           }
           int reg1 = rr1 / 1000;
           int reg2 = rr2 / 1000;
           int out1 = rr1 % 1000;
           int out2 = rr2 % 1000;
           jfile << "In random_correlate region 1 " << reg1 << " rr1 " << rr1 << " region 2 " << reg2 << " rr2 " << rr2 << endl;
           if ((reg1 == reg2) || (reg1 == out2) || (reg2 == out1)) {
              jfile << "in random_correlate neighbour detected reg 1 " << reg1 << " reg 2  " << reg2 << endl;
              continue;
           }
           if ((s1 != 0) && (s2 != 0)) { //neither edge is empty
               first = this->edgeVals[rr1];
               second = this->edgeVals[rr2];
               //printFLTVect(vectfile,first);
               //printFLTVect(vectfile,second);
               if (s1 < s2) {
                   dinterp = equalize_vector(first , second);
                   s3 = dinterp.size();
                   if (s3 == 0) {
                       cout << "Error: interpolation routine failed" << endl;
                       corr = -6;
                   }
                   else {
                       corr = correlate_Eqvector(dinterp, second);
                   }
                   if (corr <= -2) {
                       continue;
                   }
               }
               else if (s1 > s2) {
                   dinterp = equalize_vector(second, first);
                   s3 = dinterp.size();
                   if (s3 == 0) {
                       cout << "Error: interpolation routine failed" << endl;
                       corr = -6;
                   }
                   else {
                       corr = correlate_Eqvector(dinterp, first);
                   }
                   if (corr <= -2) {
                       continue;
                   }
               }
               else {
                   corr = correlate_Eqvector(first, second);
                   s3 = first.size() - second.size();
               }
               jfile <<  " r1 " << r1 << " rr1 " << rr1 <<" s1 " << s1 << " r2 " << r2 << " rr2 " << rr2 << " s2 " << s2 << " s3 " << s3 << " correlate " << corr << endl;
               kfile << corr << endl;
           } // end of if test for empty edge
        count++;
        } // end of while loop
    }

    //function to return the correlation of two vectors
    FLT correlate_Eqvector(vector<FLT> vector1, vector<FLT> vector2, bool lzero=true) {
        ofstream jfile;
        jfile.open(logpath +"/correlateEqvector.out",ios::app);
        FLT result;
        jfile << " In correlateEqvector vector 1 size " << vector1.size() << " vector 2 size " << vector2.size() << endl;
        if (vector1.size() != vector2.size()){
            jfile << "error: vectors must be same length" << endl;
            return -2;
        }
        if (lzero) {
            vector1 = this->meanzero_vector(vector1);
            vector2 = this->meanzero_vector(vector2);
        }
        FLT vector1Norm = 0;
        FLT vector2Norm = 0;
        FLT vector12Product = 0;
        unsigned int vectSize = vector1.size();
        for (unsigned int  i = 0; i < vectSize; i++) {
            vector1Norm +=  vector1[i]*vector1[i];
            vector2Norm += vector2[i]*vector2[i];
            vector12Product += vector1[i]*vector2[i];
        }
        if (vector2Norm == 0) {
            return -4;
        }
        if (vector1Norm == 0) {
            return -3;
        }
        vector1Norm = sqrt(vector1Norm);
        vector2Norm = sqrt(vector2Norm);
        result = vector12Product / (vector1Norm*vector2Norm);
        jfile << "result = " << result << endl;
        return result;

    } //end of function correlateEqvector

    vector <FLT> equalize_vector (vector<FLT> svector, vector<FLT> lvector) {
        ofstream kfile;
        kfile.open(logpath + "/eqVector.out",ios::app);
        vector <FLT> result;
        int sSize, lSize = 0;
        FLT sStep, lStep = 0;
        result.resize(0);
        sSize = svector.size();
        lSize = lvector.size();
        sStep = 1.0 / (1.0 * (sSize-1));
        lStep = 1.0 / (1.0 * (lSize-1));
        FLT start = 0;
        FLT finish = 0;
        FLT value = 0;
        FLT delta = 0.0000001;
        int marker = 0;
        for (int i=0; i<sSize-1; i++) { // walk along the short vector
            start = i*sStep;
            finish = (i+1)*sStep + delta;
            while ((marker  < lSize) && (marker*lStep < finish)) { //walk along the long vector
                    value = (svector[i+1]*(marker*lStep - start) + svector[i]*(finish - marker*lStep))/sStep;
                    result.push_back(value);
                    marker++;
            }
        }
        if (marker != lSize){
           kfile <<  " lSize " << lSize << " sSize " << sSize << " count " << marker <<  " not filled" << endl;
            result.resize(0);
            return result;
        }
        else {
             kfile << " l size " << lSize << " sSize " << sSize << "resSize " << result.size() << endl;
            // result.push_back(svector[sSize - 1]);
            return result;
        }
    } //end of function equalize_vector


    vector <FLT> resize_vector (vector<FLT> vect, int rSize) {
        ofstream qfile;
        qfile.open(logpath + "/resizeVector.out",ios::app);
        vector <FLT> result;
        int vSize = 0;
        FLT vStep, rStep = 0;
        result.resize(0);
        vSize = vect.size();
        //if the input vector is already the correct length then return the input vector
        if (vSize == rSize) {
            result = vect;
            return result;
        }

        vStep = 1.0 / (1.0 * (vSize-1));
        rStep = 1.0 / (1.0 * (rSize-1));
        FLT start = 0;
        FLT finish = 0;
        FLT value = 0;
        FLT delta = 0.0000001;
        int stepCount = 0;

        if (vSize < rSize) {
            for (int i=0; i<vSize-1; i++) { // walk along the standard size steps and interpolate
                start = i*vStep;
                finish = (i+1)*vStep + delta;
                while ((stepCount  < rSize) && (stepCount*rStep < finish + delta)) { //walk along the long
                    value = (vect[i+1]*(stepCount*rStep - start) + vect[i]*(finish - stepCount*rStep))/vStep;
                    result.push_back(value);
                    stepCount++;
                }
            }
            if (stepCount != rSize){
                qfile <<  " Error: rSize " << rSize << " sSize " << vSize << " count " << stepCount <<  " not filled" << endl;
                result.resize(0);
                return result;
            }
            else {
                qfile << " r size " << rSize << " sSize " << vSize << " result size " << result.size() << endl;
                return result;
            }
        }
        else {
            stepCount = 0; //counts increments of the standard vector
            int marker = 0; //counts increments of the variable vector
            for (int i=0; i<rSize-1; i++) { // walk along the standard steps and sum
                start = i*rStep;
                finish = (i+1)*rStep + delta;
                FLT vectCount = 0.0; //counts how many of the vector points are between standard points
                FLT value = 0.0; //stores the values of the vectors to weight the standard points
                while ((marker  < vSize) && (marker*vStep < finish + delta)) { //walk along the long
                    value += vect[i];
                    vectCount += 1.0;
                    marker++;
                }
                if (vectCount < 0.1) {
                    continue;
                }
                else {
                    result.push_back(value / vectCount);
                    stepCount++;
                }
            }
            result.push_back(vect[vSize-1]);
            stepCount++;
        }
        if (stepCount != 21){
            qfile <<  " Error: vSize " << vSize << " rSize " << rSize << " count " << stepCount <<  " not filled" << endl;
            result.resize(0);
            return result;
        }
        else {
            qfile << " vSize " << vSize << " rSize " << rSize << "result size " << result.size() << endl;
            return result;
        }
    } //end of function equalize_vector

// method to populate vector of Polygon boundary curves
    void populateBoundPolygon(bool first) {
        for (unsigned int i=0;i<NUMPOINTS;i++)
        {
            this->curvedBoundary[i] = this->polygonBoundary(i,first);
        }
    }

// method to populate vector of Polygon boundary curves
    void populateBoundCurve(bool first) {
        for (unsigned int i=0;i<NUMPOINTS;i++)
        {
            this->curvedBoundary[i] = this->roundBoundary(i,first);
        }
    }


// method to round the corners of a region
    morph::BezCurvePath<FLT> roundBoundary (int regNum, bool first)
    {
        morph::BezCurvePath<FLT> bound;
        vector<hexGeometry::point> vtxCoords;
        vector<hexGeometry::point> mtxCoords;
        int size = this->regionVertex[regNum].size();
        if (first) {
            vtxCoords = this->vCoords[regNum];
        }
        else {
            vtxCoords = this->mCoords[regNum];
        }

      //iterate over polygon vertices
        cout << " vtxCoords region " << regNum << " vtxCoords size " << vtxCoords.size() << " mCoords size " << this->mCoords[regNum].size() << " vCoords size " << this->vCoords[regNum].size() << endl;
      // code for all iterations to round corners
        for (int i=0;i<size;i++)
            {
                FLT x = (vtxCoords[i].first + vtxCoords[(i+1)%size].first)/2.0;
                FLT y = (vtxCoords[i].second + vtxCoords[(i+1)%size].second)/2.0;
                hexGeometry::point p;
                p.first = x; p.second = y;
                mtxCoords.push_back(p);
            }
        cout<< endl;
      //now create the BezCurvePaths
        for  (int i = 0; i < size;i++)
            {
                std::pair<FLT, FLT> ma, mb, va;
                ma = hGeo->point2pair(mtxCoords[((i-1)+size)%size]);
                mb = hGeo->point2pair(mtxCoords[i]);
                va = hGeo->point2pair(vtxCoords[i]);
                morph::BezCurve<FLT> bc(ma,mb,va);
                bound.addCurve(bc);
            }

        return bound;
    }//end of roundBoundary method

    /*!
     * method to determine a polygonal bezCurvePath around a region
     */
    morph::BezCurvePath<FLT> polygonBoundary (int regNum, bool first) {
        morph::BezCurvePath<FLT> bound;
        vector<hexGeometry::point> vtxCoords;
        int size = this->regionVertex[regNum].size();
        if (first) {
            vtxCoords = this->vCoords[regNum];
        }
        else {
            vtxCoords = this->mCoords[regNum];
        }

        //iterate over polygon vertices
        cout << " vtxCoords region " << regNum << " vtxCoords size " << vtxCoords.size() << " mCoords size " << this->mCoords[regNum].size() << " vCoords size " << this->vCoords[regNum].size() << endl;
          // code for all iterations to round corners
          //now create the BezCurvePaths
        for (int i = 0; i < size;i++) {
            std::pair<FLT, FLT> va, vb;
            hexGeometry::point pa = this->vCoords[regNum][((i-1)+size)%size];
            hexGeometry::point pb = this->vCoords[regNum][i];
            va = hGeo->point2pair(pa);
            vb = hGeo->point2pair(pb);
            morph::BezCurve<FLT> bc(va,vb);
            bound.addCurve(bc);
        }
        return bound;
    }//end of polygonBoundary method


    /*!
     * method to determine the line segments around a region
     */
    vector<hexGeometry::lineSegment> polygonSides (int regNum, bool first) {
        vector<hexGeometry::lineSegment> segments;
        vector<hexGeometry::point> vtxCoords;
        int size = this->regionVertex[regNum].size();
        if (first) {
            vtxCoords = this->vCoords[regNum];
        }
        else {
            vtxCoords = this->mCoords[regNum];
        }

        //iterate over polygon vertices
        cout << " vtxCoords region " << regNum << " vtxCoords size " << vtxCoords.size() << " mCoords size " << this->mCoords[regNum].size() << " vCoords size " << this->vCoords[regNum].size() << " size in polygonSides " << size << endl;
          // code for all iterations to round corners
          //now create the line segments that make the sides
        for (int i = 0; i < size;i++) {
            hexGeometry::point pa = this->vCoords[regNum][((i-1)+size)%size];
            hexGeometry::point pb = this->vCoords[regNum][i];
            cout << "pa " << pa.first << " , " << pa.second << " pb " << pb.first << " , " << pb.second << endl;
            segments.push_back(hGeo->hexGeometry::createLineSegment(pa,pb));
        }
        return segments;
    }//end of polygon side method

    bool hexInRegion(int regNum, const morph::Hex h, bool first=true) {
        bool result = true;
        hexGeometry::point p;
        p.first = h.x;
        p.second = h.y;
        vector<hexGeometry::lineSegment> segments;
        segments = this->polygonSides(regNum,first);
        vector<hexGeometry::lineSegment>::iterator ptrS;
        for (ptrS = segments.begin(); ptrS < segments.end(); ptrS++) {
            if (this->hGeo->pointLeftSegment(p, *ptrS)) {
                result = false;
                break;
            }
        }
        return result;
    }




// to find baryCentre of a region
    std::pair<FLT,FLT> baryCentre(int regNum) {
        std::pair<FLT,FLT> result;
        int bsize = regionHex[regNum].size();
        for (auto& h : regionHex[regNum]) {
             result.first += h.x;
             result.second += h.y;
        }
        result.first = result.first / (1.0 * bsize);
        result.second = result.second / (1.0 * bsize);
        return result;
    }


// returns the shortest distace from the seed point to the region boundary
// use of Creg makes it only relevant to unmorphed code
      FLT min_radius(int regNum, bool bary=true) {
          std::pair<FLT, FLT>  barycentre;
          std::pair<FLT, FLT> boundHex;
          boundHex.first = 0.0;
          boundHex.second = 0.0;
          // morphing needs to work with centres, not barycentres
          if (bary) {
              barycentre.first = this->centroids[regNum].first;
              barycentre.second = this->centroids[regNum].second;
          }
          else {
              barycentre.first =  this->centres[regNum].first;
              barycentre.second =  this->centres[regNum].second;
          }
          FLT minradius = 100000.0;
          int count = 0;
          FLT boundDist;
          for (auto h : this->regionHex[regNum]) {
              if (Creg[h.vi] > 0) {
                   count++;
                   boundHex.first = h.x,
                   boundHex.second = h.y;
                   boundDist = getdist(boundHex,barycentre);
                   if (boundDist < minradius) {
                       minradius = boundDist;
                   }
              }
          }
          cout << " minradius count of boundary hexes for region " << regNum << " is " << count << endl;
          return minradius;
      }

     FLT max_radius(int regNum, bool bary=true) {
         std::pair<FLT, FLT> boundHex;
         std::pair<FLT, FLT>  barycentre;
         boundHex.first = 0.0;
         boundHex.second = 0.0;
         //morphing needs to work with centres, not barycentre
         if (bary) {
             barycentre.first = this->centroids[regNum].first;
             barycentre.second = this->centroids[regNum].second;
         }
         else {
             barycentre.first = this->centres[regNum].first;
             barycentre.second = this->centres[regNum].second;
         }
         FLT maxradius = -100000.0;
         int count=0;
         for (auto h : this->regionHex[regNum]) {
                 count++;
                 boundHex.first = h.x,
                 boundHex.second = h.y;
                 FLT boundDist = getdist(boundHex,barycentre);

                 if (boundDist > maxradius)
                     maxradius = boundDist;
         }
         cout << " maxradius count of boundary hexes for region " << regNum << " is " << count << endl;
         return maxradius;
      }

    vector<FLT> nnRegional(int regNum) {
        vector<FLT> result;
        for (auto h : regionHex[regNum]) {
            result.push_back(nn[h.vi]);
        }
        return result;
    }

    /*
     * takes the boundary of a region and sorts it by angle
     * produces sorted vectors for
     * 1. hex indices of the hexes of the sorted boundary
     * 2. NN values of the hexes in angular order around the boundary
     * 3. The phi value for each hex in sorted order
     * all vectors have the same length and the index of the vector refers to
     * the same hex in all three arrays. It differs from the sorting in
     * renewDissect because the bounary is not dissected
     */
    void sortRegionBoundary(int regNum) {
        vector<int> regionBoundary;
        vector<int> irB;
        vector<FLT> rB;
        regionBoundary.resize(0);
        rB.resize(0);
        irB.resize(0);
        int bsize = regionBound[regNum].size();
        this->sortedBoundary[regNum].clear();
        for (auto& h : this->regionBound[regNum]) {
            FLT angle = h.phi;
            rB.push_back(angle);
            regionBoundary.push_back(h.vi);
        }
        irB = sort_indexes(rB); //indices after sort on theta
        for (int i=0; i< bsize; i++) {
            this->sortedBoundary[regNum].push_back(regionBoundary[irB[i]]);
            this->sortedBoundaryPhi[regNum].push_back(rB[irB[i]]);
            this->sortedBoundaryNN[regNum].push_back(NN[regNum][regionBoundary[irB[i]]]);
        }
    }

    /* sectorize over radius analogue version
     * unlike the old region.h we cannot directly fill from the NN field because this is different when we integrate
     * using the methods in region.h where the h.vi refer to the global positions of the hexes and are unique
     * and the use from main programs that use ksSolver for each region
     */
    vector <FLT> sectorize_reg_radius (int regNum, int _numSectors, int beginAngle, int endAngle, vector<FLT> fieldVal) {
        ofstream dfile ("logs/sectorRadius.txt",ios::app);
        int numSectors = _numSectors/2;
        vector <FLT>  radiusNN;
        vector <FLT> normalNN;
        vector <int> radiusCount;
        radiusCount.resize(numSectors,0);
        radiusNN.resize(numSectors,0);
        FLT startradius, finishradius, radiusInc; //sector radii
        FLT minradius = min_radius(regNum, true);
        FLT maxradius = max_radius(regNum, true);
        dfile << "region " << regNum << " minradius used " << minradius << " maxradius used " << maxradius <<endl;
        radiusInc = maxradius /(1.0*numSectors);
        FLT startAngle, finishAngle, angleInc; //sector angles
        angleInc = 2*PI/(1.*numSectors);
        startAngle = (beginAngle)*angleInc;
        finishAngle = (endAngle)*angleInc;
 //int size = (int) this->regionHex[regNum].size();
 // to normalise the NN field
        FLT field = 0;
        for (unsigned int i=0; i<fieldVal.size(); i++){
            normalNN.push_back(fabs(fieldVal[i]));
        }
        dfile << "after normalise field "<< field << endl;
//for (int i=0;i<size;i++)
// dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

        for (int k=0;k<numSectors;k++) {
            startradius = (k*radiusInc);
            finishradius = (k+1)*radiusInc;
            int count = 0;
            for (auto h : this->regionHex[regNum]) {
                if (h.phi >= startAngle && h.phi < finishAngle) {
                    if (h.r >= startradius && h.r <  finishradius) {
                        radiusCount[k]++;
                        radiusNN[k] += normalNN[count];
                    } //end of if on radius
                } //end of if on angleSector
                count++;
            } //end of loop over the hexes in the region
        if (radiusCount[k] == 0)
            radiusNN[k] = 0.0;
        else
            radiusNN[k] = radiusNN[k]/radiusCount[k];
        dfile << " region " << regNum << " startradius "<<startradius<<"  finishradius "<<finishradius<< " radiusNN " << radiusNN[k] << endl;

        }//end loop over regions

        dfile << endl;
        return radiusNN;

    } //end of function sectorize_radius

    /* sectorize over radius adapted for digital output
     * unlike the old region.h we cannot directly fill from the NN field because this is different when we integrate
     * using the methods in region.h where the h.vi refer to the global positions of the hexes and are unique
     * and the use from main programs that use ksSolver for each region
     */
    vector <int> sectorize_reg_Dradius (int regNum, int _numSectors, int beginAngle, int endAngle, vector<FLT> fieldVal) {
       ofstream dfile ( "logs/sectorRadius.txt",ios::app);
       int numSectors = _numSectors/2;
       vector <int>  radiusNN;
       radiusNN.resize(numSectors,0);
       vector <int> radiusCount;
       vector <FLT> radiusHold;
       vector <FLT> normalNN;
       radiusCount.resize(numSectors,0);
       radiusHold.resize(numSectors,0);
       FLT startradius, finishradius, radiusInc; //sector radii
       FLT maxradius = max_radius(regNum,true);
       FLT minradius = min_radius(regNum,true);
       dfile << "region " << regNum << " maxradius used " << maxradius << " minradius used " << minradius <<endl;
       radiusInc = maxradius /(1.0*numSectors);
       FLT startAngle, finishAngle, angleInc; //sector angles
       angleInc = 2*PI/(1.*numSectors);
       startAngle = (beginAngle)*angleInc;
       finishAngle = (endAngle)*angleInc;
        for (unsigned int i=0; i<fieldVal.size(); i++){
            normalNN.push_back(fabs(fieldVal[i]));
        }
// to normalise the NN field
       normalNN = meanzero_vector(normalNN);
       FLT epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));
 //for (int i=0;i<size;i++)
 // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;
      for (int k=0;k<numSectors;k++) {
         startradius = (k*radiusInc);
         finishradius = (k+1)*radiusInc;
         int count = 0;
         for (auto h : this->regionHex[regNum]) {
            if (h.phi >= startAngle && h.phi < finishAngle) {
                if (h.r >= startradius && h.r < finishradius) {
                    radiusCount[k]++;
//radiuscc[k] += this->cc[this->regionHex[regNum][i]];
                   radiusHold[k] += normalNN[count];
                } //end of if on radius
            } //end of if on angleSector
            count++;
         } //end of loop over hexes in and individual region
      }//end of loop over all regions

      dfile << "after creation of sectorized field region Dregion " << regNum <<  endl;

      for (int k=0;k<numSectors;k++){
         startradius = (k*radiusInc);
         finishradius = (k+1)*radiusInc;
         if (radiusCount[k] == 0) {
            radiusNN[k] = 2;
            continue;
         }
         radiusHold[k]  = radiusHold[k]  / (1.*radiusCount[k]);
         if (radiusHold[k] > epsilon)
             radiusNN[k] = 1;
         else if (radiusHold[k] < -epsilon)
             radiusNN[k] = -1;
         else
             radiusNN[k] = 0;

         dfile << " region " << regNum <<" startradius "<<startradius<<"  finishradius "<<finishradius<< " DradiusNN " << radiusNN[k] << endl;
      }//end loop over sectors
      dfile << endl;
      return radiusNN;

   } //end of function sectorize_radius



    /* function to count the hexes in sectors of a region via angular sectors
     * unlike the old region.h we cannot directly fill from the NN field because this is different when we integrate
     * using the methods in region.h where the h.vi refer to the global positions of the hexes and are unique
     * and the use from main programs that use ksSolver for each region
     */
    vector <FLT> sectorize_reg_angle (int regNum, int numSectors, int beginradius, int endradius, vector<FLT> fieldVal) {
    //std::pair<FLT,FLT> diff; //difference between seed point and CoG of region
        ofstream cfile ("logs/sectorAngle.txt",ios::app);
        vector <FLT> angleNN; //average value of cc in each sector
        vector <FLT> normalNN;
        vector <int> angleCount; //number of hexes in each sector
        angleNN.resize(numSectors,0);
        angleCount.resize(numSectors,0);
        FLT startAngle, endAngle, angleInc; //sector angles
        FLT startradius, finishradius,radiusInc;
        FLT minradius = min_radius(regNum,true);
        FLT maxradius = max_radius(regNum,true);
        radiusInc = maxradius/ (1.0*numSectors);
        startradius = beginradius*radiusInc;
        finishradius = endradius*radiusInc;
        angleInc = 2*PI/(1.*numSectors);
       cfile << "region " << regNum << " maxradius used " << maxradius << " minradius used " << minradius <<endl;
// to normalise the NN field
//int size = (int) this->regionHex[regNum].size();
        for (unsigned int i=0; i<fieldVal.size(); i++){
            normalNN.push_back(fieldVal[i]);
        }
        /*
        for (auto h : this->regionHex[regNum]){
            normalNN.push_back(fieldVal[h.vi]);
        }
        */
        normalNN = meanzero_vector(normalNN);
//for (int i=0;i<size;i++)
//   cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

        for (int k=0;k<numSectors;k++) {
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if ((k+1) == numSectors)
               endAngle = 2*PI;

            int count = 0;
            for (auto h : this->regionHex[regNum]) {
                if ( h.r  >= startradius && h.r < finishradius) {
                    if (h.phi >= startAngle && h.phi < endAngle) {
                        angleCount[k]++;
//angle[k] += this->[this->regionHex[regNum][i]];
                        angleNN[k] += normalNN[count];
                    }//end if on angle
//cfile << setw(5) << angleVal[i]  <<"  ";
                }//end if on radius
                count++;
            }//end loop over all hexes in a region
        }//end loop on all sectors
        //cfile << "after creation of sectorized field region angle " << regNum << " number of hexes " << count <<  endl;

        angleNN = meanzero_vector(angleNN);
        for (int k=0;k<numSectors;k++){ //calculate the average angle in the sector
            startAngle = k*angleInc;
            endAngle = k*angleInc;
            if (angleCount[k] != 0)
                angleNN[k] = angleNN[k]/(1.*angleCount[k]);
            else
                angleNN[k] = -999.999;
//write out values
            cfile << " region " << regNum <<" startAngle "<< startAngle << "  endAngle "<< endAngle << " angleNN " << angleNN[k] << endl;
        }//end loop on sectors

        cfile << endl;
        return angleNN;

    } //end of function sectorize_region

    /* function to count the hexes in sectors of a region via angular sectors digital version
     * unlike the old region.h we cannot directly fill from the NN field because this is different when we integrate
     * using the methods in region.h where the h.vi refer to the global positions of the hexes and are unique
     * and the use from main programs that use ksSolver for each region
     */
    vector <int> sectorize_reg_Dangle (int regNum, int numSectors, int beginradius, int endradius, vector<FLT> fieldVal) {
        ofstream cfile ("logs/sectorAngle.txt",ios::app);
 //std::pair<FLT,FLT> diff; //difference between seed point and CoG of region
        vector <int> angleNN; //digitized value of NN in each sector
        vector <FLT> angleHold;
        vector <FLT> normalNN;
        vector <int> angleCount; //number of hexes in each sector
        angleNN.resize(numSectors,0);
        angleHold.resize(numSectors,0);
        angleCount.resize(numSectors,0);
        FLT startAngle, endAngle, angleInc; //sector angles
        FLT startradius, finishradius,radiusInc;
        FLT maxradius = max_radius(regNum,true);
        FLT minradius = min_radius(regNum,true);
        radiusInc = maxradius/ (1.0*numSectors);
        startradius = beginradius*radiusInc;
        finishradius = endradius*radiusInc;
        angleInc = 2*PI/(1.*numSectors);
        cfile << "region " << regNum << " maxradius used " << maxradius << " minradius used " << minradius <<endl;
// to normalise the NN field
        for (unsigned int i=0; i<fieldVal.size(); i++){
            normalNN.push_back(fieldVal[i]);
        }
        normalNN = meanzero_vector(normalNN);
        FLT epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));

        for (int k=0;k<numSectors;k++) {
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if ((k+1) == numSectors)
            endAngle = 2*PI;
            cfile << " start of numSectors loop " << k << endl;
            int count = 0;
            for (auto &h : this->regionHex[regNum]) {
               // cfile << " start of region index  loop " << count <<  " hex " << h.vi << " h.phi " << h.phi << " h.r " << h.r << endl;

                if (h.r >= startradius && h.r < finishradius) {
                    if (h.phi >=startAngle && h.phi < endAngle) {
//angleVal.push_back(h.phi);
                    angleCount[k]++;
//angle[k] += this->[this->regionHex[regNum][i]];
                    angleHold[k] += normalNN[count];
                    }//end if on angle
//cfile << setw(5) << angleVal[angleCount[k]]  <<" region  " << regNum << endl;
                }//end if on radius
                count++;
            } //end if over hexes in a region
        }//end if on over all sectors
        cfile << "after creation of sectorized field region Dangle " << regNum <<  endl;

        angleHold = meanzero_vector(angleHold);

        for (int k=0;k<numSectors;k++){
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if (angleCount[k] == 0) {
                angleNN[k] = 2;
                continue;
            }
//angleHold[k] = angleHold[k] / (1.*angleCount[k]);
            if (angleHold[k] > epsilon)
                angleNN[k] = 1;
            else if (angleHold[k] < epsilon)
                angleNN[k] = -1;
            else
                angleNN[k] = 0;

            cfile << " region " << regNum <<" startangle  " << startAngle << "  endAngle  "<< endAngle << " DangleNN " << angleNN[k] << endl;
        } //end loop over sectors
        return angleNN;

    } //end of function sectorize_region digital version

  //method for comparing hexes by angle
  bool hexcompare(morph::Hex h1, morph::Hex h2)
  {
    bool result;
    result = (h1.phi >= h2.phi);
    return result;
  }


  //method to renew a region after rounding
  void renewRegion(int regNum, list<morph::Hex> hexen)
  {
    this->regionHex[regNum].clear();
    cout << "regionHex" << regNum << " size " << this->regionHex[regNum].size() << endl;
    for (auto& h : hexen) //repopulate regionHex
    {
      this->regionHex[regNum].push_back(h);
    }
  }

  //method to renew polars and boundary
    void renewBoundary(int regNum, list<morph::Hex> hexList, bool lNoKS=false) {
        this->regionBound[regNum].clear();
        if (lNoKS) {
            for (auto h : hexList) {
                if (this->Creg[h.vi] > 0) {
                    this->regionBound[regNum].push_back(h);
                }
            }
        }
        else {
            for (auto h : hexList) {
                if (h.boundaryHex()) {
                    this->regionBound[regNum].push_back(h);
                }
            }
        }
        cout << " region " << regNum << " bound size " <<regionBound[regNum].size() << endl;
        return;
    }

    int regionBessel(int regNum, int radialOrder, int angularOrder, FLT phase, FLT Dn=1.0) {
        int numHexes = 0;
        for (auto h : this->regionHex[regNum]) {
            FLT radius = h.r * (1.0 * Dn + 1.0);
            FLT angle = 0;
            //cout << "In regionBessel after r and before phi " << h.vi << endl;
            if (h.phi >= phase)
                angle = h.phi - phase;
            else
                angle = h.phi + 2.0 * PI - phase;
            //cout << "In regionBessel hex before Bessel call  " << h.vi << endl;
            NN[regNum].push_back(boost::math::cyl_bessel_j(radialOrder, radius) * cos(angularOrder * angle));
            //cout << "In regionBessel hex after Bessel call  " << h.vi << endl;
            numHexes++;
        }
        return numHexes;
    }

    void renewCentroids(int regNum) {
    // fill the centroids
        for (unsigned int i=0; i<NUMPOINTS; i++) {
            this->centroids[i] = this->baryCentre(i);
        }
    }

/*!
 * to clear the edges map
 */
    void edges_clear() {
        this->edges.clear();
        this->edgeIndex.clear();
    }



/*!
 * this redissects the boundary of a region
 */
    void renewDissect(int regNum, int morphNum) {
        string str = to_string(morphNum);
        string dissectFile =  this->logpath + "/renewDissect" + str + ".out";
        ofstream hhfile (dissectFile ,ios::app );
        vector<int> regionBoundary; //contains the boundary h.vi values equivalent to regionBoundary in dissectBoundary
        vector<FLT> rB; //contains the boundary thetas
        vector<int> irB; //holds the indicies of rB in angular order
        int bsize = this->regionBound[regNum].size();
        int vsize = regionVertex[regNum].size();
        hhfile << "region " << regNum << " size " << bsize << " vertex size " << vsize << endl;
        this->sortedBoundary[regNum].clear();
        vector<FLT> vertexAngle;
        // clear edges
        // fill rB with the polar angles of the boundary hexes
        for (auto& h : regionBound[regNum]) {
            FLT angle = h.phi;
            rB.push_back(angle);
            regionBoundary.push_back(h.vi);
            //cout << "region " << regNum <<" theta boundary " << angle << " boundary index " << h.vi << endl;
        }
        irB = sort_indexes(rB); //indices after sort on theta
        //hhfile << " renewDissect " << " irB size " << irB.size()  << " rB size " << rB.size() << " bsize " << bsize << endl;
        for (int i=0; i< bsize; i++) {
            this->sortedBoundary[regNum].push_back(regionBoundary[irB[i]]);
        }
        //print out the polar angles of the line segments
        for (int i=0;i<vsize;i++){
            vertexAngle.push_back(this->radialAngles[regNum][i]);
        }
        for (int i=0;i<vsize;i++){
            hhfile << "Vertex angle " << i << " = " << vertexAngle[i] << " radial angle " << this->radialAngles[regNum][i] << " ";
        }
        hhfile << endl;
        //write the indices in phi order
       for (int i = 0; i < bsize; i++)
       {
           //hhfile << " boundHex " << i << " hex " << regionBoundary[irB[i]] << " angle " << rB[irB[i]] << endl;
       }
        int offset = 0;
        int idissect = 0;
        vector<vector<int>> ihE; //Different from first dissect, we now know how many edges
        ihE.resize(vsize);
        // find the offset to the first vertex
        while ((rB[irB[offset]] < vertexAngle[0]) && (offset<bsize)) //while the angle is less than the first vertex
        {
            //hhfile << " offset inside " << offset << " vertexAngle " << vertexAngle[0] << " hex angle " << rB[irB[offset]]<< " hex " << irB[offset] << endl;
            offset++;
        }
        //hhfile << " offset outside " << offset << " hex angle " << rB[irB[offset]] << " vsize " << vsize << " hex " << irB[offset] << endl;
        idissect = offset;
        for (int i=1; i< vsize; i++)
        {
            //hhfile << "head of segment loop " << i << " idissect  " << idissect << " vertexAngle " << vertexAngle[i%vsize] << " hex " << regionBoundary[irB[idissect%bsize]] << endl;
            while ((rB[irB[idissect%bsize]] < vertexAngle[i%vsize]) && (idissect <= bsize))
            {
          //      hhfile << " filling edge loop " << i-1 << " index " << idissect << " angle " << rB[irB[idissect%bsize]] << " hex " << regionBoundary[irB[idissect%bsize]] << endl;
                ihE[i-1].push_back(regionBoundary[irB[idissect%bsize]]);
                idissect++;
            }
        }
        //hhfile << " just before  vsize end loop angle " << vertexAngle[0] + 2*PI << endl;
        while ((rB[irB[idissect%bsize]] < vertexAngle[0] + 2*PI) && (idissect< bsize+offset-1))
        {
          //  hhfile << " filling end edge loop 2 " << idissect << " angle " << rB[irB[idissect%bsize]] << " hex " << regionBoundary[irB[idissect%bsize]] <<  endl;
            ihE[vsize-1].push_back(regionBoundary[irB[idissect%bsize]]);
            idissect++;
        }

        for (int iregion = 0;iregion<vsize;iregion++)
        {
            int edgeOuter = this->regionList[regNum][iregion];
            //hhfile << "edgeOuter " << edgeOuter << " region " << regNum << endl;
            if (edgeOuter > -1)
            {
                std::pair<int,int> keypair(regNum,edgeOuter);
                int keyint = this->pair2int(keypair,this->base);
                std::pair <int, vector<int>> p1(keyint,ihE[iregion]);
            //    hhfile << "region " << regNum << " outer " << iregion << " ihE " << endl;
                printIntVect(dissectFile, ihE[iregion]);
                if (this->edges.insert(p1).second) {
                    this->edgeIndex.push_back(keyint);
                }
                else {
                    cout << "duplicate entry in region" << regNum << " outer " << edgeOuter << " keyp " << keyint << endl;
                }
            }
        } // end filling edges loop

        //now print out the edges
        /*
        vector<int>::iterator ptr;
        for (ptr = regionList[regNum].begin(); ptr < regionList[regNum].end(); ptr++)
        {
            std::pair <int,int> klook(regNum, *ptr);
            int k = pair2int(klook,this->base);
            if (edges.count(k) == 0) {
                continue;
            }
            else {
                int sizeij = edges[k].size();
                 hhfile << regNum << " " << *ptr << " sizeij " << sizeij << " key " << k << endl;
                 printIntVect(dissectFile,edges[k]);
            }
        }*/
        hhfile << "Printing out edges " << endl;
        for (auto itr = edges.begin(); itr != edges.end(); itr++) {
            hhfile << " key " << itr->first << " list " << endl;
            printIntVect(dissectFile, itr->second);
        }
        hhfile << " renewDissect size of edges " << this->edges.size() << " size of edgeIndex " << this->edgeIndex.size() << endl;
    } //end of method renewDissect


    // function to renew correlate matching edges
    FLT renewcorrelate_edges(int regNum,  const int morphNum, bool lZero=true)
    {
        FLT result = 0;
        string str = to_string(morphNum);
        ofstream corrfile(this->logpath + "/correlate" + str + ".data",ios::app);
        ofstream edgefile(this->logpath + "/edgeCorrelations" + str + ".txt",ios::app);
        edgefile << " in morphed edge correlation routine "<<endl;
        vector<FLT> first;
        vector<FLT> second;
        vector<FLT> dinterp;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        FLT NNmean1, NNmean2;
        NNmean1 = this->meanNN(regNum);
        edgefile << " mean of NN in region " << regNum << "is " << NNmean1 << endl;
        for (auto j = this->regionList[regNum].begin(); j < this->regionList[regNum].end(); j++)  {
            if (*j == -1) continue;
            NNmean2 = this->meanNN(*j); //find the mean of NN in the region
            edgefile << " j iteration " << *j << " region " << regNum <<" NNmean2 " << NNmean2 << endl;
            first.resize(0);
            second.resize(0);
            dinterp.resize(0);
            std::pair<int,int> edgePair1(regNum,*j);
            int edgeIndex1 = this->pair2int(edgePair1,this->base);
            std::pair<int,int>  edgePair2(*j,regNum);
            int edgeIndex2 = this->pair2int(edgePair2,this->base);
            edgefile << "region " << regNum << " outer " << *j << " index1 " << edgeIndex1 << " index2 " << edgeIndex2 << endl;
            int count1 = 0;
            int count2 = 0;
            FLT correlationValue = -3;
            FLT ratio;
            //create two vectors from the opposing edges
            for (auto itr = this->edges[edgeIndex1].begin(); itr != this->edges[edgeIndex1].end();itr++)
            {
                first.push_back(this->NN[regNum][*itr]);
                count1++;
            }
            first = this->meanzero_vector(first, NNmean1);
            for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++)
            {
                second.push_back(this->NN[*j][*itr]);
                count2++;
			}
            second = this->meanzero_vector(second, NNmean2);
            //edges on either side are traversed in opposite directions
            std::reverse(second.begin(),second.end());
            for (unsigned int  i = 0; i < first.size(); i++) {
                if (isnan(first[i]))
                     edgefile << "Nan at " << i << " in first" << endl;
            }
            for (unsigned int  i = 0; i < second.size(); i++) {
                if (isnan(second [i]))
                     edgefile << "Nan at " << i << " in second" << endl;
            }
            //do this if the vectors are the same size
            if (first.size() == second.size() && second.size()*first.size() != 0)
            {
                correlationValue = this->correlate_Eqvector(first, second, lZero);
                ratio = 1.0;
                result += fabs(correlationValue);
                countResult++;
                corrfile <<  correlationValue << "  " <<  endl;
                edgefile << " if 1 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            } //end of code if both edges are equal
            //do this if the first vector is bigger than the second
            else if (first.size() > second.size() && first.size()*second.size() != 0)
            {
                dinterp = this->equalize_vector(second,first);
                correlationValue = this->correlate_Eqvector(dinterp, first, lZero);
                ratio = 1.0 * second.size() / (1.0 * first.size());
                if (correlationValue > -2) {
                    result += fabs(correlationValue);
                    countResult++;
                    corrfile <<  correlationValue << "  " <<  endl;
                }
                edgefile << " if 2 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            }
            //do this if the first vector is smaller than the second
            else if (first.size() < second.size() && first.size()*second.size() != 0)
            {
                dinterp = this->equalize_vector(first,second);
                correlationValue = this->correlate_Eqvector(dinterp, second, lZero);
                ratio = 1.0 * first.size() / (1.0 * second.size());
                if (correlationValue > -2) {
                    result += fabs(correlationValue);
                    countResult++;
                    corrfile <<  correlationValue << "  " <<  endl;
                }
                edgefile << " if 3 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            }
            //something has gone wrong!
            else {
                edgefile  << " if 4 region " << regNum  <<" count1 " << first.size() << " count2 " << second.size() << " j " << *j << " dinterp is zero  "   << endl;
            }
            vector<FLT> tempVect = this->resize_vector(first,21);
            std::pair<int,vector<FLT>> p1(edgeIndex1,tempVect);
            edgeNN.insert(p1);

        } //end of loop over the edges in the region
        edgefile << " countResult "<<countResult<< " entries in edgeIndex " << this->edgeIndex.size() << " entries in edges " << this->edges.size() << endl;
        edgefile << endl;
        result = result / (countResult * 1.0);
        edgefile.close();
        corrfile.close();
        return result;
    } //end of function renewcorrelate_edges

    // function to renew correlate matching edges for special functions
    FLT Srenewcorrelate_edges(int regNum,  const int morphNum, bool lZero=true)
    {
        FLT result = 0;
        string str = to_string(morphNum);
        ofstream corrfile(this->logpath + "/correlate" + str + ".data",ios::app);
        ofstream edgefile(this->logpath + "/edgeCorrelations" + str + ".txt",ios::app);
        edgefile << " in morphed edge correlation routine "<<endl;
        vector<FLT> first;
        vector<FLT> second;
        vector<FLT> dinterp;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        FLT NNmean1, NNmean2;
        NNmean1 = this->meanNN(regNum);
        //NNmean1 = 0;
        edgefile << " mean of NN in region " << regNum << "is " << NNmean1 << endl;
        for (auto j = this->regionList[regNum].begin(); j < this->regionList[regNum].end(); j++)  {
            if (*j == -1) continue;
            NNmean2 = this->meanNN(*j); //find the mean of NN in the region
            //NNmean2 = 0;
            edgefile << " j iteration " << *j << " NNmean2 " << NNmean2 << endl;
            first.resize(0);
            second.resize(0);
            dinterp.resize(0);
            std::pair<int,int> edgePair1(regNum,*j);
            int edgeIndex1 = this->pair2int(edgePair1,this->base);
            std::pair<int,int>  edgePair2(*j,regNum);
            int edgeIndex2 = this->pair2int(edgePair2,this->base);
            edgefile << "region " << regNum << " outer " << *j << " index1 " << edgeIndex1 << " index2 " << edgeIndex2 << endl;
            int count1 = 0;
            int count2 = 0;
            FLT correlationValue = -3;
            FLT ratio;
            int edge1Size = edges[edgeIndex1].size();
            for (int i=0; i<edge1Size; i++)
            {
                if (isnan(this->NN[regNum][i]))
                     edgefile << "Nan  edge hex 1 " << i <<  endl;
                first.push_back(this->NN[regNum][i]);
                count1++;
            }
            first = this->meanzero_vector(first, NNmean1);
            int edge2Size = edges[edgeIndex2].size();
            for (int i=0; i<edge2Size; i++)
            {
                if (isnan(this->NN[*j][i]))
                     edgefile << "Nan edge 2 edge hex " << i << endl;
                second.push_back(this->NN[*j][i]);
                count2++;
            }
            second = this->meanzero_vector(second, NNmean2);


            std::reverse(second.begin(),second.end());
            for (unsigned int  i = 0; i < first.size(); i++) {
                if (isnan(first[i]))
                     edgefile << "Nan at " << i << " in first" << endl;
            }
            for (unsigned int  i = 0; i < second.size(); i++) {
                if (isnan(second [i]))
                     edgefile << "Nan at " << i << " in second" << endl;
            }
            if (first.size() == second.size() && second.size()*first.size() != 0)
            {
                correlationValue = this->correlate_Eqvector(first, second, lZero);
                ratio = 1.0;
                result += fabs(correlationValue);
                countResult++;
                corrfile <<  correlationValue << "  " <<  endl;
                edgefile << " if 1 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            } //end of code if both edges are equal
            else if (first.size() > second.size() && first.size()*second.size() != 0)
            {
                dinterp = this->equalize_vector(second,first);
                correlationValue = this->correlate_Eqvector(dinterp, first, lZero);
                ratio = 1.0 * second.size() / (1.0 * first.size());
                result += fabs(correlationValue);
                countResult++;
                corrfile <<  correlationValue << "  " <<  endl;
                edgefile << " if 2 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            }
            else if (first.size() < second.size() && first.size()*second.size() != 0)
            {
                dinterp = this->equalize_vector(first,second);
                correlationValue = this->correlate_Eqvector(dinterp, second, lZero);
                ratio = 1.0 * first.size() / (1.0 * second.size());
                if (correlationValue > -2) {
                    result += fabs(correlationValue);
                    countResult++;
                    corrfile <<  correlationValue << "  " <<  endl;
                }
                edgefile << " if 3 region " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << " ratio " << ratio << endl;
            }
            else {
                edgefile  << " if 4 region " << regNum  <<" count1 " << first.size() << " count2 " << second.size() << " j " << *j << " dinterp is zero  "   << endl;
                //corrfile << " -1.5 " << ratio << endl;
            }
        } //end of loop over the edges in the region
        edgefile << " countResult "<<countResult<< " entries in edgeIndex " << this->edgeIndex.size() << " entries in edges " << this->edges.size() << endl;
        edgefile << endl;
        result = result / (countResult * 1.0);
        edgefile.close();
        corrfile.close();
        return result;
    } //end of function renewcorrelate_edge
}; // DRegion
