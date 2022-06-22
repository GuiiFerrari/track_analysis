{
TFile* file1 = TFile::Open("suma_a_in.root");
TFile* file2 = TFile::Open("suma_a_out.root");


TH1F* histodata1 = (TH1F*) file1->Get("/ThCM_zprof/ThetaCM_hist_120");
TH1F* histodata2 = (TH1F*) file2->Get("/ThCM_zprof/ThetaCM_hist_120");
//TH1F* histodata3 = (TH1F*) file1->Get("AngD_Ex_2.25_MeV");
//TH1F* histodata4 = (TH1F*) file1->Get("AngD_Ex_20.25_MeV");


	int nbins = histodata1->GetNbinsX();
	float Solid[90] ;
	float Ruth[90];

	for(int bin =1; bin<nbins-1; bin++){

		Double_t binCenter = histodata1->GetXaxis()->GetBinCenter(bin);
		//Double_t binCenter_up = histodata1->GetXaxis()->GetBinCenter(bin);
		Double_t dOmega = 2.0*3.141592*( cos((bin-1)*2*3.141592/180.) - cos((bin)*2*3.141592/180.) );
		Solid[bin] = 	dOmega;
		//Double_t binEY = histodata1->GetBinError(bin);
		
		//cout<<binCenter<<"  "<<binY<<"  "<<binEY<<endl;

		}
	
    	Solid[89] = Solid[88];
	double m1 = 4;
	double m2 = 17;
	double z1 = 2;
	double z2 = 9;
	double K1 = 5.61;


	for(int bin =1; bin<nbins; bin++){

		Double_t binCenter = histodata1->GetXaxis()->GetBinCenter(bin);
		Double_t binY1 = histodata1->GetBinContent(bin);
		Double_t binEY1 = histodata1->GetBinError(bin);
		
		Double_t binY2 = histodata2->GetBinContent(bin);
		Double_t binEY2 = histodata2->GetBinError(bin);
		double classic = 10.*pow(z1*z2*197./137.,2.0)/( pow(4*K1*m2/(m1+m2),2.0)* pow(sin(0.5*binCenter*3.14159/180),4.0) );
		
		cout<<binCenter<<"  "<<binY1/Solid[bin]<<"  "<<binEY1/Solid[bin]<<"  "<<binY2/Solid[bin]<<"  "<<binEY2/Solid[bin]<<"  "<<classic<<endl;

		}


}

