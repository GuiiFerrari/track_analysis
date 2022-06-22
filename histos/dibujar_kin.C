{











TFile f("17f_aa_kin.root",  "RECREATE");

ifstream fER("17f_a_30.txt");//from LISE++,
		
		Double_t l1=0, l2=0, l3=0, l4=0, l5=0;
		Int_t model=1;
		
		Double_t X1[9999];
		Double_t Y1[9999];
		Double_t Z1[9999];
		int i = 0;
		for (string line; getline(fER, line);) {
			stringstream parse_die(line);
			vector<Double_t> iRE;
			parse_die >> l1 >> l2 >> l3;
			X1[i]=l1; Y1[i]=l2; Z1[i]=l3;
			i++;	
		}
		fER.close();
		
		
		
		TGraph *gr1 = new TGraph(9999,Y1,Z1);

	gr1->Draw("AC");

ifstream fER2("17f_a_25.txt");//from LISE++,
		
		 l1=0, l2=0, l3=0, l4=0, l5=0;
		
		
		Double_t X2[9999];
		Double_t Y2[9999];
		Double_t Z2[9999];
		i = 0;
		for (string line; getline(fER2, line);) {
			stringstream parse_die(line);
			
			parse_die >> l1 >> l2 >> l3;
			X2[i]=l1; Y2[i]=l2; Z2[i]=l3;
			i++;	
		}
		fER2.close();
		
		
		
		TGraph *gr2 = new TGraph(9999,Y2,Z2);

	gr2->Draw("AC");

ifstream fER3("17f_a_14.txt");//from LISE++,
		
		 l1=0, l2=0, l3=0, l4=0, l5=0;
		
		
		Double_t X3[9999];
		Double_t Y3[9999];
		Double_t Z3[9999];
		i = 0;
		for (string line; getline(fER3, line);) {
			stringstream parse_die(line);
			
			parse_die >> l1 >> l2 >> l3;
			X3[i]=l1; Y3[i]=l2; Z3[i]=l3;
			i++;	
		}
		fER3.close();
		
		
		
		TGraph *gr3 = new TGraph(9999,Y3,Z3);

	gr3->Draw("AC");





	f.WriteObject(gr1,"17f_kin_30");
	f.WriteObject(gr2,"17f_kin_25");
	f.WriteObject(gr3,"17f_kin_14");



//gr2->Draw("AC");
//f->WriteObject(gr2,"20ne_p1");
//gr3->Draw("AC");
//f->WriteObject(gr3,"20ne_pd");
//gr4->Draw("AC");
//f->WriteObject(gr4,"56ni_p1_400");
//mg->Add(gr2);
//mg->Draw();
//gr1->Wite();

//a->Draw("DSSD2:-atan(pos_x_dssd2/353)*180/3.1415+32.5>>kin+line"," ener_dssd2>0 ", "colz");
//gDirectory->Append(gr1);





//gr1->Wite();




//gr2->Wite();
//gr3->Wite();
//gr4->Wite();

f.Close(); //cerramos el archivo de salida






}
