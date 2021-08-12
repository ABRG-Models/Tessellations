/*
  ksSolver class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * creates a hexGrid given a boundary curve and solves
 * the KS equations
 *
 */
#include <morph/tools.h>
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
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
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

class ksSolver
{
public:
// list of objects visible to member functions
    int scale;
	FLT xspan;
    int n;
    FLT ds;
	FLT nnInitialOffset = 1.0;
	FLT ccInitialOffset = 2.5;
	FLT boundaryFalloffDist = 0.024;
    pair<FLT,FLT> seedPoint;
	BezCurvePath<float> bound;
	string logpath;
    vector<vector<int> > N; // hex neighbourhood
    vector<FLT> NN, CC; //hold the field values for each he
    morph::HexGrid* Hgrid;
// empty constructor
    ksSolver(){};
// constructor with HexGrid passed in
    ksSolver (morph::HexGrid*  Hgrid, std::string logpath) {
        this->Hgrid = Hgrid;
        this->seedPoint = this->Hgrid->computeCentroid(this->Hgrid->hexen);
        this->logpath = logpath;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        this->seedPoint = this->Hgrid->computeCentroid(this->Hgrid->hexen);
        this->ds = this->Hgrid->hexen.begin()->d;
        n = 0;
        n = this->Hgrid->num();
        afile << " max x " << this->Hgrid->getXmax(0.0) << " min x " << this->Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << this->Hgrid->num() << endl;
        afile << "after creating HexGri ds = "<< this->ds << endl;
        cout << "after creating HexGri ds = "<< this->ds << endl;
        afile << "seed point.x " << seedPoint.first << " seed point.y " << seedPoint.second << endl;
        n = this->Hgrid->num();
        afile << " max x " << this->Hgrid->getXmax(0.0) << " min x " << this->Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
      // check the order numbering in hexen
        this->N.resize(n);
        this->setHexType();
        this->NN.resize(n);
        CC.resize(n);
        afile << "after alloc NN and CC" <<endl;
        pair<FLT, FLT> centroid = set_kS_polars(this->seedPoint);
        cout << " end of ksSolver from file " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
        cout << " end of ksSolver from file " << " centroid x " << centroid.first << " centroid y " << centroid.second << endl;
    }; // end of ksSolver constructor

// constructor with radius passed in for solving on radial boundaries
    ksSolver (int scale, FLT xspan, string logpath, float radius, pair<float, float> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        afile << "after creating HexGrid with radius "<< radius << endl;
        Hgrid->setCircularBoundary(radius, seedPoint, false);
        afile << "after setting boundary on  H " << Hgrid->num() << " centroid.x " << Hgrid->boundaryCentroid.first << " centroid.y " << Hgrid->boundaryCentroid.second << endl;
        afile << "after setting boundary on  H " << Hgrid->num() << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in circular constructor" << " n = " << n <<endl;
        N.resize(n);
        for (auto &h : this->Hgrid->hexen){
            this->N[h.vi].resize(6);
            if (!HAS_nE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][0] = h.vi;
            }
            else {
                this->N[h.vi][0] = Hgrid->d_ne[h.vi];
            }

            if (!HAS_nNE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][1] = h.vi;
            }
            else {
                this->N[h.vi][1] = Hgrid->d_nne[h.vi];
            }

            if (!HAS_nNW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][2] = h.vi;
            }
            else {
                this->N[h.vi][2] = Hgrid->d_nnw[h.vi];
            }

            if (!HAS_nW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][3] = h.vi;
            }
            else {
                this->N[h.vi][3] = Hgrid->d_nw[h.vi];
            }

            if (!HAS_nSW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][4] = h.vi;
            }
            else {
                this->N[h.vi][4] = Hgrid->d_nsw[h.vi];
            }

            if (!HAS_nSE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][5] = h.vi;
            }
            else {
                this->N[h.vi][5] = Hgrid->d_nse[h.vi];
            }
        } //end of loop over HexGri this->setHexType();
        /*
        Hgrid->computeDistanceToBoundary();
        for (auto &h : Hgrid->hexen) {
            afile  << "dist to bdry " << h.distToBoundary << " for " << h.vi << endl;
        }
        */
        this->NN.resize(n);
        CC.resize(n);
        afile << "after alloc NN and CC" <<endl;
        pair<FLT, FLT> centroid = set_kS_polars(this->seedPoint);
        cout << " end of ksSolver circle radius " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
        cout << " end of ksSolver circle radius " << " centroid x " << centroid.first << " centroid y " << centroid.second << endl;
    }; // end of ksSolver constructor

// Constructor with boundary passed in
    ksSolver (int scale, FLT xspan, string logpath, BezCurvePath<float> bound, pair<FLT,FLT> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        Hgrid->setBoundary(bound,false);
        reverse_y();
        afile << "after setting boundary on  H " << Hgrid->num() << " centroid.x " << Hgrid->boundaryCentroid.first << " centroid.y " << Hgrid->boundaryCentroid.second << endl;
        afile << "seed point.x " << seedPoint.first << " seed point.y " << seedPoint.second << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
      // check the order numbering in hexen
        N.resize(n);
        this->setHexType();
        this->NN.resize(n);
        this->CC.resize(n);
        afile << "after alloc NN and CC" <<endl;
        pair<FLT, FLT> centroid = set_kS_polars(this->seedPoint);
        cout << " end of ksSolver bezCurvePath " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
        cout << " end of ksSolver bezCurvePath " << " centroid x " << centroid.first << " centroid y " << centroid.second << endl;
    }; // end of ksSolver constructor


    void setHexType() {
        for (auto &h : this->Hgrid->hexen){
            this->N[h.vi].resize(6);
            if (!HAS_nE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][0] = h.vi;
            }
            else {
                this->N[h.vi][0] = Hgrid->d_ne[h.vi];
            }

            if (!HAS_nNE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][1] = h.vi;
            }
            else {
                this->N[h.vi][1] = Hgrid->d_nne[h.vi];
            }

            if (!HAS_nNW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][2] = h.vi;
            }
            else {
                this->N[h.vi][2] = Hgrid->d_nnw[h.vi];
            }

            if (!HAS_nW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][3] = h.vi;
            }
            else {
                this->N[h.vi][3] = Hgrid->d_nw[h.vi];
            }

            if (!HAS_nSW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][4] = h.vi;
            }
            else {
                this->N[h.vi][4] = Hgrid->d_nsw[h.vi];
            }

            if (!HAS_nSE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][5] = h.vi;
            }
            else {
                this->N[h.vi][5] = Hgrid->d_nse[h.vi];
            }
        } //end of loop over HexGri this->setHexType();
        /*
        Hgrid->computeDistanceToBoundary();
        for (auto &h : Hgrid->hexen) {
            afile  << "dist to bdry " << h.distToBoundary << " for " << h.vi << endl;
        }
        */
    } // end of setHexType

// Constructor with boundary passed in

// method to calculate the Laplacian
    vector<FLT> getLaplacian(vector<FLT> Q, FLT dx) {
        FLT overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<FLT> L(n,0.);
        for(auto &h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overds;
        }
        return L;
    }

		vector<FLT> chemoTaxis(vector<FLT> Q, vector<FLT> P, FLT dx) {
		vector<FLT> cT(n,0.);
          FLT overds = 1./(1.5*29.0*29.0*dx*dx);

        for (auto &h : Hgrid->hexen) {
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

	  //G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds) + NN[i]*lapC[i];

        } //matches for on i
		return cT;
	} //end of function chemoTaxis

  // function to compute the derivative
     void compute_dNNdt(vector<FLT>& inN, vector<FLT>& dNdt, FLT Dn, FLT Dchi) {
        vector<FLT> lapN(this->n,0);
		vector<FLT> cTaxis(this->n,0);
		//cout << "in compute_dNN just before laplacian" << endl;
		lapN = getLaplacian(inN,this->ds);
        cTaxis = chemoTaxis(inN,this->CC,this->ds);
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
            N2 = this->NN[h]*this->NN[h];
            dCdt[h] =  beta*N2/(1.+N2) - mu*inC[h] + Dc*lapC[h];
        }
    }

  //function to timestep coupled equations solely b.c. on the flux
    void step(FLT dt, FLT Dn, FLT Dchi, FLT Dc)
	{
        dt = dt * 2.5 / Dn;

        // cout  << "value of NN[5] start Runge " << this->NN[5] << endl;


        // Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
		{
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
            if(h.boundaryHex())
			{
                for(int j=0;j<6;j++)
				{
		            int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
                    if(N[h.vi][j] == i) {
       	                this->NN[N[h.vi][j]] = this->NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	                    this->CC[N[h.vi][j]] = this->CC[h.vi];
	                }
	            }
	         }
          }

        // 2. Do integration of NN
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
            this->compute_dNNdt (this->NN, dNdt, Dn, Dchi);
			 //cout << "in ksSolver after compute_dNNdt" << endl;
            for (int h=0; h< this->n; ++h) {
                K1[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h< this->n; ++h) {
                K2[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dNNdt (Ntst, dNdt, Dn, Dchi);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dNdt[h] * dt;
                Ntst[h] = this->NN[h] + K3[h];
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
                this->NN[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
				//this->NN[i] = i * 1.0;
            }
        }

        // 3. Do integration of B
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
            this->compute_dCCdt (this->CC, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K1[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K1[h] * 0.5 ;
            }

            /*
             * Stage 2
             */
		    this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K2[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K2[h] * 0.5;
            }

            /*
             * Stage 3
             */
            this->compute_dCCdt (Ctst, dCdt, Dc);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dCdt[h] * dt;
                Ctst[h] = this->CC[h] + K3[h];
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
                this->CC[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
            }
        }
        //cout  << "value of NN[5] end Runge " << this->NN[5] <<  " number of hexes " << this->n << endl;
    }//end step

  //function to timestep coupled equations option to set boundary to constant value
    void step(FLT dt, FLT Dn, FLT Dchi, FLT Dc, int steps, int numAdjust)
	{
        dt = dt * 2.5 / Dn;


      if ((steps%numAdjust == 0) && (steps/numAdjust != 0))
	  {
	     cout << "in numAdjust if step " << steps << endl;
	     for (auto &h : this->Hgrid->hexen)
	     {
		     //cout << "dist to bdry" << h.distToBoundary << endl;
	         if (h.distToBoundary > -0.5)
		     { // It's possible that distToBoundary is set to -1.0
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- this->boundaryFalloffDist)) );
				//cout << "bSig " << bSig << " for hex " << h.vi << endl;
                this->NN[h.vi] = (this->NN[h.vi] - this->nnInitialOffset) * bSig + this->nnInitialOffset;
                this->CC[h.vi] = (this->CC[h.vi] - this->ccInitialOffset) * bSig + this->ccInitialOffset;
		     } //end of if on boundary distance
	     }//end of loop over hexGrid
	  } //end of code applied to keep boundary conditions static

        // Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
		{
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
            if(h.boundaryHex())
			{
                for(int j=0;j<6;j++)
				{
		            int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
                    if(N[h.vi][j] == i)
					{
       	                this->NN[N[h.vi][j]] = NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	                    this->CC[N[h.vi][j]] = this->CC[h.vi];
	                }
	            }
	         }
          }
          void step(FLT dt, FLT Dn, FLT Dchi, FLT Dc);
    }//end step

    void reverse_y ()
	{
	  for (auto &h : this->Hgrid->hexen)
	    {
	      //cout << " in reverse_y " << h.vi << endl;
	      //int index = h.vi;
	      //cout << " in y reversing loop " << h.vi << endl;
	      float temp = h.y;
	     // cout << " after getting y " << temp << endl;
		  if (temp != 0.0)
		  {
	         // cout << " in y reversing loop " << endl;
	          h.y = -temp;
              this->Hgrid->d_y[h.vi] = -temp;
		  }
		}
	}

// function to give r and theta relative to region centre
    pair <FLT,FLT> set_kS_polars(pair<FLT,FLT> centre){
        pair <FLT, FLT> result;
		result.first = 0.0;
		result.second = 0.0;
        FLT xav=0;
        FLT yav = 0;
        int hexcount = 0;
        FLT maxPhi = -10.0;
        FLT minPhi = 10.0;
	    cout <<"in set polars ksSolver xcentre" << centre.first << " y_centre  " << centre.second <<endl;
        for (auto &h : this->Hgrid->hexen) {
            hexcount++;
            xav += h.x;
            yav += h.y;
        }
		//cout << "after set centres " << endl;
		if (hexcount != 0) {
        xav = xav / (hexcount*1.0);
        yav = yav / (hexcount*1.0);
		}
		else {
		  //cout << " in set_polars no hexes in region "<<endl;
		  }
//go over the region and put the hexes into bins then average
        for (auto&  h : this->Hgrid->hexen) {
			FLT angle = 0;
            FLT dx = h.x;
            FLT dy = h.y;
            h.r = sqrt((dx - centre.first)*(dx - centre.first)
			+ (dy - yav)*(dy - yav));
            if (dy >= centre.second) {
               angle =   atan2((dy - centre.second), (dx - centre.first));
			  h.phi = angle;
			  }
            else {
                angle =  2*PI + atan2((dy - centre.second), (dx - centre.first));
			    h.phi = angle;
			}
            if (angle < minPhi) {
                minPhi = angle;
            }
            if (angle > maxPhi) {
                maxPhi = angle;
            }
        }
        cout << " set_kS_polars max phi " << maxPhi << " minPhi " << minPhi << endl;
        result.first = xav - centre.first ; // barycentre
        result.second = yav - centre.second;
		cout << "centre x "<< centre.first << " centre y " << centre.second << " centreMove.x " << result.first << " centreMove.y " << result.second <<endl;
        return result;
    } //end of function set_polars

}; //end of class KSsolver
