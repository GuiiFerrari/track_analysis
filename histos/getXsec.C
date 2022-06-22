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







void getXsec() {

TGaxis::SetMaxDigits(3);  ///notacion cientifica para los ejes a partir de 10^n

	TFile* file1 = TFile::Open("suma_pnew3.root");



	for(int i=0;i<40;i++){
		ostringstream filenamedwba ;
		//filenamedwba << std::fixed;
   		//filenamedwba.precision(2);
		int w = i*10;	
		//filenamedwba<<"/Th_zprof/Theta_hist_"<<w;
        	filenamedwba<<"/Th_zprof/Theta_hist_C_"<<w;
		TH1F* histodata1 = (TH1F*) file1->Get(filenamedwba.str().data());        

		ostringstream fileener;
        	fileener<<"/Th_zprof/EBeam_hist_"<<w;
		TH1F* histodata2 = (TH1F*) file1->Get(fileener.str().data());        



		int nbins = histodata1->GetNbinsX();
		float Solid[45] ;


		for(int bin =1; bin<nbins-1; bin++){

			//Double_t binCenter = histodata1->GetXaxis()->GetBinCenter(bin);
			//Double_t binCenter_up = histodata1->GetXaxis()->GetBinCenter(bin);
			Double_t dOmega = 2.0*3.141592*( cos((bin-1)*2*3.141592/180.) - cos((bin)*2*3.141592/180.) );
			Solid[bin] = 	dOmega;
			//Double_t binEY = histodata1->GetBinError(bin);
		
			//cout<<binCenter<<"  "<<binY<<"  "<<binEY<<endl;

		}
	
    		Solid[44] = Solid[43];

		cout<<endl;
		cout<<endl;
		//cout<<'\t'<<"  "<<w<<'\t'<<histodata2->GetMean()<<"  "<<histodata2->GetRMS()<<endl;
		//cout<<endl;

		for(int bin =1; bin<nbins; bin++){

			Double_t binCenter = histodata1->GetXaxis()->GetBinCenter(bin);
			Double_t binY1 = histodata1->GetBinContent(bin);
			Double_t binEY1 = histodata1->GetBinError(bin);
		
			
			cout<<w<<"  "<<histodata2->GetMean()<<"  "<<histodata2->GetRMS()<<"  "<<binCenter<<"  "<<binY1<<"  "<<binEY1<<"  "<<Solid[bin]<<endl;

		}



        }  






}

