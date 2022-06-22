#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"


#include<string>
#include <cstdio>
#include <sstream>
#include <iostream>
#include<fstream>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include <stdlib.h> 
#include <ctime>


#define NF 90
#define NT 179
#define pi 3.1415926535
#define Nang 7
#define NangT 178


//#include "get_strength.h"

using namespace std;







void fitXsec() {

TGaxis::SetMaxDigits(3);  ///notacion cientifica para los ejes a partir de 10^n



	ostringstream filename ;
	filename <<"alpha_17F_angdis_220.dat";
        
  	ifstream  entradaXsec;
	double angle0, angle1, angle_av, ex0,  ex1, ex_av, xsec, dxsec ;
        Double_t Angle0[NF], Angle1[NF], Angle_av[NF], Ex0[NF], Ex1[NF], Ex_av[NF], Xsec[NF], dXsec[NF];
	
        entradaXsec. open(filename.str().data());
      	if(entradaXsec.fail() ){
                       cerr << "error abriendo "<< filename.str().data() << endl;
 			exit(1);
                      }  

	                
		int goodpoints = 0;
         	for(int k=0;k<NF;k++){
        		entradaXsec >> ex0 >> xsec >> dxsec >> angle0 >> angle1   ;
			//Ex0[k]=ex0; Ex1[k]=ex1;   Ex_av[k]=ex_av;
			 //Angle1[k]=ex1;  Angle_av[k]=ex_av;
			//Xsec[k]=xsec*2.0/0.4; dXsec[k]=dxsec*2.0/0.4; //because the rebinning was not taken into account before
			//cross section table is in units of mb/(sr 400keV)
			if(ex0>=75 && ex0<=103){Angle0[goodpoints]=ex0; Xsec[goodpoints]=angle0; dXsec[goodpoints]=angle1; 
			cout<<  Angle0[goodpoints] <<"  "<< Xsec[goodpoints] <<"  "<< dXsec[goodpoints]<<endl ;
			goodpoints++;
			}

			

		

	        }

                entradaXsec. close();


	



		
		TGraphErrors* Data = new TGraphErrors(goodpoints,Angle0,Xsec,0,dXsec );
		//TGraphErrors* Data = new TGraphErrors(Nang-1,ANG,XS,0,XSE );
		//TGraphErrors* Data = new TGraphErrors(Nang-2,ANG,XS,0,XSE );




		//------------------------Leemos DWBA


		ostringstream filenamedwba ;
		//filenamedwba << std::fixed;
   		//filenamedwba.precision(2);
        	filenamedwba<<"spp_220.dat";        
                
        
        
  		ifstream  entradadwba;
		double angle, c2L0,  c3L1, c4L2, c5L3 ;
        	Double_t X[NangT], YL0[NangT], YL1[NangT], YL2[NangT], YL3[NangT];
	
        	entradadwba. open(filenamedwba.str().data());
      		if(entradadwba.fail() ){
                       cerr << "error abriendo "<< filenamedwba.str().data() << endl;
 			exit(1);
                      }  

                
         	for(int k=0;k<NangT;k++){
        		//entradadwba >> angle >> c2L0 >> c3L1 >> c4L2 >> c5L3   ;
			//X[k]=angle; YL0[k]=c2L0;   YL1[k]=c3L1; YL2[k]=c4L2; YL3[k]=c5L3;						
			entradadwba >> angle >> c2L0 >> c3L1   ;
			X[k]=angle; YL0[k]=c2L0;   YL1[k]=c3L1; 
			//cout<<X[k]<<"  "<<YL0[k]<<endl;

	        }

                entradadwba.close();

                TGraph* GL0 = new TGraph(NangT,X,YL0);
                //GL0->SetTitle("L = 0");
                /*TGraph* GL1 = new TGraph(NangT,X,YL1);
                GL1->SetTitle("L = 1");
                TGraph* GL2 = new TGraph(NangT,X,YL2);
		GL2->SetTitle("L = 2");
                TGraph* GL3 = new TGraph(NangT,X,YL3);
		GL3->SetTitle("L = 3");
                TGraph* GL23 = new TGraph(NangT,X,YL2);
                GL23->SetTitle("L #geq 2");
                TGraph* GLTot = new TGraph(NangT,X,YL0);
                GLTot->SetTitle("Total");
		*/
                //GL0->Draw("AC*");
                TF1 * f0 = new TF1("f0",[&](double*x, double *par){ return par[0]*GL0->Eval(x[0]) ; }, 20, 120, 1);
                //TF1 * f1 = new TF1("f1",[&](double*x, double *par){ return par[0]*GL1->Eval(x[0]) ; }, 0, 180, 1);
                //TF1 * f2 = new TF1("f2",[&](double*x, double *par){ return par[0]*GL2->Eval(x[0]) ; }, 0, 180, 1);
                //TF1 * f3 = new TF1("f3",[&](double*x, double *par){ return par[0]*GL3->Eval(x[0]) ; }, 0, 180, 1);

                //TF1 * function = new TF1("Total",[&](double*x, double *par){ return par[0]*GL0->Eval(x[0]); }, 20, 120, 1);
	//TF1 * function = new TF1("Total",[&](double*x, double *par){ return par[0]*GL0->Eval(x[0]) + par[1]*GL1->Eval(x[0]) + par[2]*GL2->Eval(x[0]); }, 0, 10, 3);

       
        	//f0->SetLineColor(kRed); GL0->SetLineColor(kRed); GL0->SetLineWidth(2); GL0->SetFillColor(kWhite);
        	//f1->SetLineColor(kBlue); GL1->SetLineColor(kBlue); GL1->SetLineWidth(2); GL1->SetFillColor(kWhite);
		//f2->SetLineColor(kGreen); GL2->SetLineColor(kGreen); GL2->SetLineWidth(2); GL2->SetFillColor(kWhite);
        	//f2->SetLineColor(kGreen); GL23->SetLineColor(kGreen); GL23->SetLineWidth(2); GL23->SetFillColor(kWhite);
        	//f3->SetLineColor(kOrange+1); GL3->SetLineColor(kOrange+1); GL3->SetLineWidth(2); GL3->SetFillColor(kWhite);
        	//function->SetLineColor(kMagenta); GLTot->SetLineColor(kMagenta); GLTot->SetLineWidth(2); GLTot->SetFillColor(kWhite);


        	//f0->SetParLimits(0, 0, 50);
        	//function->SetParLimits(1, 0, 50);
        	//function->SetParLimits(2, 0, 50);
        	//function->SetParLimits(3, 0, 50);

		//f0->FixParameter(0, 1);
		//function->FixParameter(3, 0);
		//function2->FixParameter(3, 0);
	
		//-----------------------------


		//--------------Hacemos el Fit

		//Data->Fit(function2,"MEQ0"); //solo para calcular errores
        	Data->Fit(f0,"ME"); //here is the magic!
        	Data->SetTitle("Exp. data");
        	Data->SetMarkerStyle(8);
		Data->SetFillColor(kWhite);

        	//Data->Draw("AP");

		//f0->Draw("L");

        	Double_t K0 = f0->GetParameter(0);
        	cout<<K0<<endl;




}

