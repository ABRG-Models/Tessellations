/*!
 * If COMPILE_PLOTTING is defined at compile time, then include the display and
 * plotting code. I usually put all the plotting code inside #defines like this so
 * that I can compile a version of the binary without plotting, for parameter searches
 * in which I am only going to be saving out HDF5 data.
 */
#ifdef COMPILE_PLOTTING
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/ColourMap.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#endif

#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/Config.h>
#include <morph/Scale.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>

//! Helper function to save PNG images with a suitable name
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, morph::Visual& v)
{
    std::stringstream ff1;
    ff1 << logpath << "/" << name<< "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}

using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using namespace morph;
using namespace std;

#ifdef COMPILE_PLOTTING
class vtxVisual : public morph::VisualModel
{
public:
    // class data
    std::vector<std::vector<morph::Vector<FLT,3>>> vtx;
    unsigned int size;
    FLT linewidth = 0.00390625;
    FLT radiusFixed = 0.00390625;
    VBOint idx = 0;
    morph::Vector<FLT,3> uz = {1.0f, 1.0f, 1.0f};
    // class constructors
    vtxVisual(GLuint sp, GLuint tsp, std::vector<std::vector<morph::Vector<FLT, 3>>> & _vtx, FLT line, FLT radius)
    {
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->vtx = _vtx;
        this->size = _vtx.size();
        this->linewidth = line;
        this->radiusFixed = radius;
    }
    //class methods
    void initializeVertices() {
        this->initv();
        cout << " in initializeVertices after call to initv " << endl;
    }

    // Draw vertices for the net's actual locations
    void initv() {
    // Discs at the net vertices
        cout << "entering intitv " << endl;
        morph::Vector<FLT,3> puckthick = { 0, 0, 0.002 };
        morph::Vector<FLT,3> linethick = {0, 0, 0.001};
        morph::Vector<FLT,3> colStart = {0.0, 0.0, 0.0};
        morph::Vector<FLT,3> colEnd = {0.0, 0.0, 0.0};
        for (unsigned int j=0; j < this->size; ++j) {
            unsigned int jsize = vtx[j].size();
            for (unsigned int i=0; i < jsize; i++) {
                cout << " initv vtx " << " i " << i << " j " << j << " is " << vtx[j][i] << endl;
                this->computeTube (this->idx,
                      this->vtx[j][i]+puckthick,
                      this->vtx[j][i]-puckthick,
                      morph::Vector<FLT,3>({1,0,0}), morph::Vector<FLT,3>({0,1,0}),
                      colStart, colEnd,
                      radiusFixed, 16);
            }
            cout << "after inner loop for computeTube " << endl;
        }
        cout << "after connect tube" << endl;
        // Connection lines
        for (unsigned int j=0; j < size; j++) {
            unsigned int jsize = vtx[j].size();
            for (unsigned int i=0; i < jsize; i++) {
                morph::Vector<FLT, 3> c1 = this->vtx[j][i] + linethick;
                morph::Vector<FLT, 3> c2 = this->vtx[j][(i+1)%jsize] + linethick;
                this->computeLine (idx, c1, c2, this->uz, colStart, colEnd, linewidth, linewidth/4);
            }
        }
    }
}; //end of class vtxVisual
#endif

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << " supply the json file as arg[1]" << endl;
      return -1;
    }
    string jsonfile = argv[1];
    //  open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
#ifdef SINGLE
        float dt = conf.getFloat("dt",0.0001);
        float Dn = conf.getFloat("Dn",1.0);
        float Dchi = conf.getFloat("Dchi",0.0);
        float Dc = conf.getFloat("Dc",0.3);
        float xspan = conf.getFloat("xspan",5.0);
        float boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float nnInitialOffset = conf.getFloat("nnInitialOffet", 1.0);
        float ccInitialOffset = conf.getFloat("ccInitialOffset",2.5);
#else
        double dt = conf.getDouble("dt",0.0001);
        double Dn = conf.getDouble("Dn",1.0);
        double Dchi = conf.getDouble("Dchi",0.0);
        double Dc = conf.getDouble("Dc",0.3);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffet", 1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
#endif
    int numSectors = conf.getInt("numsectors",12);
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numprint = conf.getInt("numprint",100);
    string logpath = conf.getString("logpath", "./logsMorph") ;
    string iter = conf.getString("iter","0");
    bool LfixedSeed = conf.getBool("LfixedSeed",0);
    bool LDn = conf.getBool("LDn",0);
    //bool overwrite_logs = conf.getBool("overwrite_logs",true);
    bool skipMorph  = conf.getBool("skipMorph",false);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    unsigned int numpoints = conf.getInt("numpoints",41);
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
    // adjust the number of steps according to the Dn number
    numsteps = numsteps * floor(sqrt(36.0/Dn));
    numprint  = numprint * floor(sqrt(36.0/Dn));
    // adjust the time step for the Dn values
    dt = dt * sqrt(Dn/36.0);
    //set up a vtxVisual pointer
    vtxVisual* cv;

    unsigned int seed;
    if (LfixedSeed) {
        seed = 1;
    }
    else {
        seed = time(NULL);
    }

    // A ra2yyndo2yym uniform generator returning real/FLTing point types
    morph::RandUniform<FLT> ruf(seed);
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");
    ofstream degfile1 (logpath + "/degree1.data");
    ofstream degfile2 (logpath + "/degree2.data");
    ofstream degfile3 (logpath + "/degree3.data");


// initialise DRegion class setting scale
    DRegion M(scale,xspan,logpath,numpoints); //create tessellation
    M.setCreg(); //set counts to identify inner boundaries
    M.setInternalBoundary(); //set internal boundaries
    cout << "before dissect_boundary " << endl;
    vector<std::pair<FLT,FLT>> cGravity;
    cGravity = M.dissectBoundary(); //dissect region boundary
    cout << "Edges size = " << M.edges.size() << endl;
	M.setRadialSegments(); //set the radial segments for regions
    int inReg = 0;
    inReg = M.setInnerRegion(); //set mask array for inner regions
    int rC = 0;
    //check for correct number of inner regions
    for (unsigned int j=0; j<numpoints;j++) {
        if (!M.innerRegion[j]) rC++;
    }
    if (inReg != rC) {
        cout << "Error: setInnerRegion returns " << inReg << " but outerRegions " << rC << endl;
        return -1;
    }
    else {
        cout << "Success: setInnerRegion no of outer regions " << inReg << endl;
    }
    cout << "after first setRadialSegments " << endl;
// include the analysis methods
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph0 "<< j <<endl;
    }
    cout << "first setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "just before populating the inner boundary Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
        //M.renewCentroids(j);
    }
    cout << "just before calculating regionSize " << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
	cout << "just after populating the regions from the ksSolver vector"<<endl;
    //clear global edges map
    M.edges_clear();
    cout << "just before renewDissect first time" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,0);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;
 #ifdef COMPILE_PLOTTING
// now draw the intial tesselation
    FLT hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
    // Parameters from the config that apply only to plotting:
    const unsigned int plotevery = conf.getUInt ("plotevery", 10);
    // Should the plots be saved as png images?
    const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", false);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);
    cout << "just before new Visual object" << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual * v1;
    v1 = new morph::Visual(win_width, win_height, "Tessellation0 ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1->backgroundWhite();
    //orthographic projection
    v1->ptype = morph::perspective_type::orthographic;
    v1->ortho_bl = {-1.0f, -1.0f};
    v1->ortho_tr = {1.0f, 1.0f};
    // You can tweak the near and far clipping planes
    v1->zNear = 0.001;
    v1->zFar = 20;
    // And the field of view of the visual scene.
    v1->fov = 45;
    // You can lock movement of the scene
    v1->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1->setZDefault (conf.getFloat ("z_default", 0.0f));
    v1->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1->scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    //make it current
    v1->setCurrent();

    // if using plotting, then set up the render clock
    steady_clock::time_point lastrender = steady_clock::now();
// to draw the tessellation
// first convert M.vCoords into a morph::Vect
    cout << "just beofere creating vtxVector" << endl;
    std::vector<std::vector<morph::Vector<FLT,3>>> vtxVector;
    vtxVector.resize(M.vCoords.size());
    cout << "just after creating vtxVector" << endl;
    for (unsigned int j=0;j<numpoints;j++) {
        for (unsigned idx = 0; idx<M.vCoords[j].size(); idx++) {
            vtxVector[j].push_back(M.hGeo->point2vect3(M.vCoords[j][idx]));
        }
    }
    cout << "after filling vtxVector " << endl;
    // now instantiate vtxVisual
    cv = new vtxVisual(v1->shaderprog, v1->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 0", {0.0f, 1.1f, 0.0f});
    //v1->addVisualModel(cv);
    cout << "before rendering v1 " << endl;
    v1->render();
    cout << "after rendering v1 " << endl;
    std::stringstream frame;
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    //frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation0", 0, *v1);


#endif

// initialise the fields
    string fname = logpath + "/first.h5";
    cout << "just before first data read morph 0 "<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++)
        {
		    std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
	  }
      else
	  {
        for (unsigned int j=0;j<numpoints;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        FLT choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

 #ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    Vector<float> spatOff;
    float xzero = 0.0f;

    // A. Offset in x direction to the left.
    xzero -= 0.5*0.008;
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    Scale<FLT,float> cscale; cscale.do_autoscale = true;
    unsigned int Agrid[numpoints];

    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
        if (M.innerRegion[j]) {
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Agrid[j] = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                                    v1->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
        }//end of loop on inner regions
    }//end of loop over regions
#endif

    // begin morph0 time stepping loop
    for (int i=0;i<numsteps;i++) {
   	    for (unsigned int j = 0;j<numpoints;j++) { //loop over all regions, only step internal ones
            if (M.innerRegion[j]) {
                S[j].step(dt, Dn, Dchi, Dc);
            }
        }
#ifdef COMPILE_PLOTTING
        if ((i % plotevery) == 0) {
            for (unsigned int j=0; j<numpoints;j++) {
                if (M.innerRegion[j]) { //only display inner regions
                    vector<FLT> regionNN;
                    regionNN = L.normalise(S[j].NN);
                    VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v1->getVisualModel (Agrid[j]);
                    avm->updateData (&regionNN);
                    avm->clearAutoscaleColour();
                    if (saveplots) {
                        if (vidframes) {
                            savePngs (logpath, "nn0", framecount, *v1);
                             ++framecount;
                        }
                        else {
                             savePngs (logpath, "nn0", i , *v1);
                        }
                    }
                }
            }
        }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 17000) { // 17 is about 60 Hz
            glfwPollEvents();
            v1->render();
            lastrender = steady_clock::now();
        }
    }
  //  cerr << "Ctrl-c or press x in graphics window to exit.\n";
  //  v1->keepOpen();
#endif
    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before first data write morph 0" << endl;
    morph::HdfData fdata(fname);
	for (unsigned int j=0;j<numpoints;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        fdata.add_contained_vals(ccst,S[j].CC);
        fdata.add_contained_vals(nst,S[j].NN);
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;
    gfile << endl << "analysis on first morphing iteration " << endl;

    //declaration of variables needed for analysis
    vector <int> radiusDVector;
    vector <int> angleDVector;
	vector <FLT> angleVector;
    vector <FLT> radiusVector;
    int degreeRadius;
    int degreeAngle;
	FLT tempArea = 0.0;
	FLT tempPerimeter = 0.0;
    int angleOffset = 0;
    int radiusOffset = 0;
	FLT avDegreeRadius = 0.0;
	FLT avDegreeAngle = 0;
	FLT occupancy = 0.0;
    int countRegions = 0;
    FLT avAbsCorrelation = 0.0;
    const int max_comp = numpoints*3;
    for (unsigned int j=0;j<numpoints;j++) {
        M.NN[j] = S[j].NN; //write the NN fields to the DRegion array
        if (M.innerRegion[j]){
			int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);

            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M.meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			    radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			    newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			    if (newdegreeRadius > degreeRadius) {
                    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  cout << "just after renewcorrelate_edges morph1 " << endl;
          //avAbsCorrelation = M.correlate_edges(0);
		  M.random_correlate(max_comp, 1);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (unsigned int j=0;j<numpoints;j++) {
	        if (M.innerRegion[j]) {
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              avAbsCorrelation += M.renewcorrelate_edges(j,1);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
      if (countRegions == 0) {
          cout << "Error zero regionss counted in second analysis morph 0" << endl;
          return -1;
      }
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
        cout << "Regions counted in first analysis loop " << countRegions << endl;
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;

        // write the edge vectors all interpolated to a uniform size
        std::map<int, vector<FLT>>::iterator ptr;
        std::string sideNN = logpath + "/edgeNN.data";
        int ecount = 0;
        for (ptr = M.edgeNN.begin(); ptr != M.edgeNN.end(); ptr++) {
            vector<FLT> tempVect;
            tempVect = M.normaliseVect(ptr->second);
            M.printFLTVect(sideNN, tempVect);
            cout << "edgeNN key " << ptr->first << " element " << ecount << " vector size " << ptr->second.size() << " edges size " << M.edges[ptr->first].size() <<  endl;
            ecount++;
        }
//end of integration after solving on polygonal regions via ksSolver

    /*
     * now create an array of morphed regions, first morph
     */

    if (skipMorph) return 0;
    cout << "just before setting curved boundaries first morph" <<endl;
// section for solving on the curved boundaries
    S.resize(0);
    //set radius for creating circular regions also calculate the adjusted Dn values
    vector<FLT> DnVal;
    DnVal.resize(numpoints,0.0);
    vector<FLT> DchiVal;
    DchiVal.resize(numpoints,0.0);
    vector<FLT> DcVal;
    DcVal.resize(numpoints,0.0);
    vector<FLT> morph0Area;
    morph0Area.resize(numpoints,0.0);
	for (unsigned int j = 0;j<numpoints;j++) {
        morph0Area[j] = M.hexArea*M.regArea(j);
    }
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(1);
    cout << "just after setting curved boundaries morph 1 " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph1   "<< j <<endl;
    }
    cout << "first morph  setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions morph 1" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "just before populating the inner boundary morph 1" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
    }
    cout << "second setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
	cout << "just after populating the ksVector"<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
	//	for (auto h : S[j].Hgrid->hexen) {
	//	    cout << "hexNumber in region " << j <<  " h.vi " << h.vi << endl;
    //    }
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] = Dchi * sqrt(morph0Area[j] / area);
            DnVal[j] = Dn *  sqrt(morph0Area[j] / area);
            DcVal[j] = Dc * sqrt(morph0Area[j] / area);
        }
        else {
            DchiVal[j] = Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }

        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue) {
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++) {
		    std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
    }
    else {
        for (unsigned int j=0;j<numpoints;j++) {
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        FLT choice = ruf.get();
		        if (choice > 0.5) {
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else {
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
	    	    } //end of if on boundary distance
	        }//end of loop over region
	    }//end of loop over all regions
    } //end of else on Lcontinue
    cout <<  "just after field creation first morph" << endl;
#ifdef COMPILE_PLOTTING
    morph::Visual * v2;
    //morph::Vector<FLT,3> offset ={1.0,0.0,0.0};
    v2 = new morph::Visual(win_width, win_height, "Tessellation1 ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v2->backgroundWhite();
    v2->ptype = morph::perspective_type::orthographic;
    v2->ortho_bl = {-1.0f, -1.0f};
    v2->ortho_tr = {1.0f, 1.0f};
    // You can tweak the near and far clipping planes
    v2->zNear = 0.001;
    v2->zFar = 20;
    // And the field of view of the visual scene.
    v2->fov = 45;
    // You can lock movement of the scene
    v2->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v2->setZDefault (conf.getFloat ("z_default", -5.0f));
    v2->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v2->scenetrans_stepsize = 0.5;
    // now instantiate vtxVisual
    cv = new vtxVisual(v2->shaderprog, v2->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 1", {0.0f, 1.1f, 0.0f});
    v2->addVisualModel(cv);
    cout << "before rendering v2 " << endl;
    v2->render();
    cout << "after rendering v2 " << endl;
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    //frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation1", 0, *v2);
    v2->setCurrent();
    cout << "before rendering v2 " << endl;
    unsigned int Bgrid[numpoints];
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
        if (M.innerRegion[j]) {
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Bgrid[j] = v2->addVisualModel (new HexGridVisual<FLT> (v2->shaderprog,
                                                                    v2->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
        }//end of loop on inner regions
    }//end of loop over region
#endif
    // begin third time stepping loop after first morph
    for (int i=0;i<numsteps;i++) {
   	for (unsigned int j = 0;j<numpoints;j++){ //loop over regions, only set internal ones
            S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
        }
#ifdef COMPILE_PLOTTING
        if ((i % plotevery) == 0) {
            for (unsigned int j=0; j<numpoints;j++) {
                if (M.innerRegion[j]) { //only display inner regions
                    vector<FLT> regionNN;
                    regionNN = L.normalise(S[j].NN);
                    VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v2->getVisualModel (Bgrid[j]);
                    avm->updateData (&regionNN);
                    avm->clearAutoscaleColour();
                }
            }
            v2->render();
            if (saveplots) {
                if (vidframes) {
                    savePngs (logpath, "nn1", framecount, *v2);
                    ++framecount;
                }
                else {
                    savePngs (logpath, "nn1", i , *v2);
                }
            }
        }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v2.render().
        //steady_clock::duration sincerender = steady_clock::now() - lastrender;
        //if (duration_cast<milliseconds>(sincerender).count() > 17000) { // 17 is about 60 Hz
        //    glfwPollEvents();
        //    v2->render();
        //    lastrender = steady_clock::now();
        //}
    }
#endif
    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before second data read " << endl;
    morph::HdfData gdata(gname);
	for (unsigned int j=0;j<numpoints;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);
     cout << " just after writing data "  << endl;
     // write the NN and CC vals for each region
      gfile << endl << "analysis on first morphing iteration " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
      for (unsigned int j=0;j<numpoints;j++) {
          M.NN[j] = S[j].NN; // first fill the DRegion NN values
              if (M.innerRegion[j]){
			      int regionCount = 0;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M.regArea(j);
                  tempPerimeter = M.renewRegPerimeter(j);

                  // digital version
                  angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M.meanzero_vector(radiusVector);

                  degreeRadius = L.find_zeroDRadius(radiusDVector);
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
				  {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				      degreeRadius = newdegreeRadius;
		          }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                  regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT


    M.edges_clear();
    // swap the radialAngles to the mCoords
    // M.swapRadialSegments(false);
    // redissect the boundaries
    cout << "just before renewDissect second time" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,1);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  M.random_correlate(max_comp,2);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (unsigned int j=0;j<numpoints;j++) {
	        if (M.innerRegion[j])
			{
	          countRegions++;
		      occupancy += M.regNNfrac(j);
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              avAbsCorrelation += M.renewcorrelate_edges(j,2);
              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on NUMPOINTs
      if (countRegions == 0) {
          cout << "Error zero regionss counted in third analysis morph 1" << endl;
          return -1;
      }
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after first morphing
//begin second morphing
    if (skipMorph) return 0;
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(0);
    S.resize(0);
    cout << "just after setting curved boundaries morph 2 " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++)
    {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the ksVector morph2 "<< j <<endl;
    }
    cout << "just after populating the ksVector"<<endl;
	// repopulate the regions
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewRegion(j,S[j].Hgrid->hexen);
        }
    }
    cout << "after repopulate regions morph2" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        }
        //M.renewCentroids(j);
    }
    cout << "after repopulate boundary morph2" << endl;
//    swap the radialAngles to the mCoords
    M.edges_clear();
    M.swapRadialSegments(false);
    // redissect the boundaries
    for (unsigned int j=0;j<numpoints;j++)
    {
        if (M.innerRegion[j]) {
            M.renewDissect(j,2);
        }
    }
    cout << "Edges size " << M.edges.size() << endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] =  Dchi * sqrt(morph0Area[j] /  area);
            DnVal[j] = Dn * sqrt(morph0Area[j] / area);
            DcVal[j] = Dc * sqrt(morph0Area[j] / area);
        }
        else {
            DchiVal[j] =  Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }
        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// now draw the intial tesselation
// initialise the fields
        cout<< "just before third data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
        string hname = logpath + "/third.h5";
        if (Lcontinue) {
             morph::HdfData hinput (hname,1);
             for (unsigned int j=0;j<numpoints;j++){
                std::string nstr = "n" + to_string(j);
                char * nst = new char[nstr.length()+1];
                std::strcpy(nst,nstr.c_str());
                std::string ccstr = "c" + to_string(j);
                char * ccst = new char[ccstr.length()+1];
                cout << "labels "<< nstr <<" , " << ccstr<<endl;
                std::strcpy(ccst,ccstr.c_str());
                hinput.read_contained_vals(nst,S[j].NN);
                hinput.read_contained_vals(ccst,S[j].CC);
                cout<< "just after input of NN and CC1"<< endl;
       //   input.close();
            }
	 }
     else {
        for (unsigned int j=0;j<numpoints;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5)
                {
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
                }
                else {
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
                }

                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
            }//end of loop over region
        }//end of loop over all regions
     } //end of else on Lcontinue
#ifdef COMPILE_PLOTTING
    morph::Visual * v3;
    v3 = new morph::Visual(win_width, win_height, "Tessellation2");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v3->backgroundWhite();
    //orthographic projection
    v3->ptype = morph::perspective_type::orthographic;
    v3->ortho_bl = {-1.0f, -1.0f};
    v3->ortho_tr = {1.0f, 1.0f};
    // You can tweak the near and far clipping planes
    v3->zNear = 0.001;
    v3->zFar = 20;
    // And the field of view of the visual scene.
    v3->fov = 45;
    // You can lock movement of the scene
    v3->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v3->setZDefault (conf.getFloat ("z_default", -5.0f));
    v3->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v3->scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    // now instantiate vtxVisual
    cv = new vtxVisual(v3->shaderprog, v3->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 2", {0.0f, 1.1f, 0.0f});
    v3->addVisualModel(cv);
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    //frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation2", 0, *v3);
    v3->setCurrent();
    //v3->render();
    unsigned int Cgrid[numpoints];

    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
        if (M.innerRegion[j]) {
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Cgrid[j] = v3->addVisualModel (new HexGridVisual<FLT> (v3->shaderprog,
                                                                    v3->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
        }//end of loop on inner regions
    }//end of loop over region
    #endif

     //begin fourth time stepping loop after second morph
    int loopsteps = 0;
    for (int i=0;i<numsteps;i++) {
        for (unsigned int j = 0;j<numpoints;j++) {
            S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
        }
#ifdef COMPILE_PLOTTING
        if((i % plotevery) == 0) {
            for (unsigned int j=0; j<numpoints;j++) {
                if (M.innerRegion[j]) { //only display inner regions
                    vector<FLT> regionNN;
                    regionNN = L.normalise(S[j].NN);
                    VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v3->getVisualModel (Cgrid[j]);
                    avm->updateData (&regionNN);
                    avm->clearAutoscaleColour();             if (saveplots) {
                        if (vidframes) {
                            savePngs (logpath, "nn2", framecount, *v3);
                             ++framecount;
                        }
                        else {
                             savePngs (logpath, "nn2", i , *v3);
                        }
                    }
                }
            }
        }
        /*
   	        // rendering the graphics. After each simulation step, check if enough time
            // has elapsed for it to be necessary to call v3.render().
            cout << " before rendering the graphics " << endl;
            v3->render();
            cout << " after rendering the graphics " << endl;
            savePngs (logpath, "nnField2", 0, *v3);
         */

        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v3.render().
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 17000) { // 17 is about 60 Hz
            glfwPollEvents();
            v3->render();
            lastrender = steady_clock::now();
        }
    }
#endif
    //
	cout << "after third time stepping loop " << loopsteps << endl;
       //nndisp.closeDisplay();
    //code run at end of timestepping
    cout << "just before the fourth data write" << endl;
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
	cout <<"after creating hdata "<< hname << endl;
	for (unsigned int j=0;j<numpoints;j++) {
		std::string nstr = "n" + to_string(j);
	   	char * nst = new char[nstr.length()+1];
		std::strcpy(nst,nstr.c_str());
		std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
		std::strcpy(ccst,ccstr.c_str());
		cout <<"in third hdf5 write "<<nstr <<" , "<<ccstr<<endl;
        hdata.add_contained_vals(ccst,S[j].CC);
        hdata.add_contained_vals(nst,S[j].NN);
	}
    //data.add_contained_vals("X",M.X[0]);
    //data.add_contained_vals("Y",M.X[1]);
    // hdata.add_val ("/Dchi", Dchi);
    // hdata.add_val ("/Dn", Dn);
    // hdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on second morphing iteration " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (unsigned int j=0;j<numpoints;j++) {
        M.NN[j] = S[j].NN;
	    if (M.innerRegion[j]) {
	        int regionCount = 0;
            gfile<<"in the degree loop" << endl;
	        //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

           // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl;
            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
			for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius) {
				    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT

//computing average values over all regions
    avDegreeAngle = 0;
	avDegreeRadius = 0;
	occupancy = 0;
    avAbsCorrelation = 0;
	tempPerimeter = 0;
	tempArea = 0;
	countRegions = 0;
	M.random_correlate(max_comp, 3);
	cout << "just after random correlate_edges morph2 " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (unsigned int j=0;j<numpoints;j++) {
	    if (M.innerRegion[j]){
            countRegions++;
            avAbsCorrelation += M.renewcorrelate_edges(j,3);
            occupancy += M.regNNfrac(j);
            cout << "after renNNfrac " << endl;
            tempArea = M.regArea(j);
            cout << " after regArea " << endl;
            tempPerimeter = M.regPerimeter(j);
            cout << " after regPerimeter " << endl;
            cout << "before sectorixe Dangle morph2 " << endl;
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            cout << "after sectorixe Dangle morph2 " << endl;
            degreeAngle = L.find_zeroDAngle(angleDVector);
            avDegreeAngle += degreeAngle;
            cout << "after sectorixe Dangle morph2 " << endl;
            //radial degree
            degreeRadius = 0;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            avDegreeRadius += degreeRadius;
            degfile3 << degreeAngle << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
         } //end of if on non-zero regions
    } //end of loop on NUMPOINTs
    if (countRegions == 0) {
        cout << "Error zero regionss counted in fourth analysis morph 2" << endl;
        return -1;
    }
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	occupancy = occupancy / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation  <<endl;
    cout << "end of program reached successfully!" << endl;
    return 0;
};
