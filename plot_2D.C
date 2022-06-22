{
   
   TFile* file1 = TFile::Open("Histoscattering_events_17F_219.root");
    //TFile* file1 = TFile::Open("ProtoOut.root");
   //TFile* file1 = TFile::Open("ATTPCOut.root");
     //TFile* file1 = TFile::Open("ACTAROut_ana.root");

 
   int evento = 4403;

   TGraph2D *d2t = (TGraph2D*) file1->Get(Form("H3D_%d",evento));
   d2t->Draw("P0");

   TGraph2D *dvertex = (TGraph2D*) file1->Get(Form("vertex_%d",evento));
   dvertex->SetMarkerColor(kBlue);
   dvertex->SetMarkerStyle(20);
   dvertex->Draw("P same");

   TGraph2D *dlpoint = (TGraph2D*) file1->Get(Form("lpoint_%d",evento));
   dlpoint->SetMarkerColor(kRed);
   dlpoint->SetMarkerStyle(20);
   dlpoint->Draw("P same");
  	
  TVectorD* numberL = (TVectorD*) file1->Get(Form("Nlines_%d",evento));
  cout<<numberL->Norm1()<<endl;

   Int_t nlines =  (Int_t) (numberL->Norm1());
   //const Int_t value = nlines;
   TPolyLine3D *line[4];
	/*line[0] = (TPolyLine3D*) file1->Get(Form("Line_%d_%d",evento,0));
	line[1] = (TPolyLine3D*) file1->Get(Form("Line_%d_%d",evento,1));
	line[0]->Draw("same");
	//line[1]->Draw("same");
	*/
  /*line[3] = (TPolyLine3D*) file1->Get(Form("Line_%d_%d",evento,3));
  line[3]->SetLineColor(kGreen);
  line[3]->Draw("same");
	*/
  for(int w=0; w<nlines; w++){
	line[w] = (TPolyLine3D*) file1->Get(Form("Line_%d_%d",evento,w));
	if(w==0){ line[w]->SetLineColor(kRed);  }
  	if(w==1){ line[w]->SetLineColor(kBlack);}
	if(w==2){ line[w]->SetLineColor(kBlue); }
	if(w==3){ line[w]->SetLineColor(kGreen); }
	line[w]->Draw("same");
	
	}

	//dvertex->Draw("P same");

 //TPolyLine3D *l1 = (TPolyLine3D*) file1->Get("Line1_1");   
   //l1->Draw(" same ");	
   /*TPolyLine3D *l2 = (TPolyLine3D*) file1->Get("Line2_132");
   l2->SetLineColor(kRed);
   l2->Draw(" same ");
 TPolyLine3D *l3 = (TPolyLine3D*) file1->Get("Line3_132");   
   l3->Draw(" same ");	
 */
 /*TPolyLine3D *l4 = (TPolyLine3D*) file1->Get("Line4_7");   
   l4->Draw(" same ");	
   */

  
   /*double xmin[2] = {0};
   double ymin[2] = {0};
   double zmin[2] = {0};
   dline->GetPoint(0, xmin[0], ymin[0], zmin[0]); 
   */
   //int ix,iy,iz; dline->GetHistogram()->GetMinimumBin(ix,iy,iz); double xmin = dline->GetHistogram()->GetXaxis()->GetBinCenter(ix); double ymin = dline->GetHistogram()->GetYaxis()->GetBinCenter(iy);   
 
   
	

   /*TGraph2D *dline2 = (TGraph2D*) file1->Get("Line2_2");
   dline2->SetLineColor(2);
			dline2->SetMarkerColor(2);
			//dline2->SetMarkerStyle(8);
			//dline2->SetMarkerSize(1);
			dline2->SetLineWidth(2);
   //dline2->Draw("same line");
*/
  
	/*
    int n = 500;
    double dt = 100.;
   // draw original line
    double xh = -7.53362 ; double xm = -7.22205;
    double yh = 62.5937 ; double ym =-0.0455406 ;
    double zh = 7.36715; double zm = 7.28445;
    TPolyLine3D *l1 = new TPolyLine3D(2*n);
    for (int i = -n; i <n;++i) {
       double t = -10 + dt*i/n;
	
       double x,y,z;
	x = xm + (xh -xm)*t;
	y = ym + (yh -ym)*t;
	z = zm + (zh -zm)*t;
	//cout<<i<<" "<<t<<"  "<<x<<"  "<<y<<"  "<<z<<endl;
       l1->SetPoint(i,x,y,z);
    }
    l1->SetLineColor(kRed);
    l1->SetLineWidth(2);
    l1->Draw("same");

   
    // draw original line
     xh = -41.3294;  xm = 43.1857;
     yh = 161.714;  ym =43.8685;
     zh = -10.8303 ;  zm = 33.9445 ;
    TPolyLine3D *l2 = new TPolyLine3D(2*n);
    for (int i = -n; i <n;++i) {
       double t = -10 + dt*i/n;
	
       double x,y,z;
	x = xm + (xh -xm)*t;
	y = ym + (yh -ym)*t;
	z = zm + (zh -zm)*t;
	//cout<<i<<" "<<t<<"  "<<x<<"  "<<y<<"  "<<z<<endl;
       l2->SetPoint(i,x,y,z);
    }
    l2->SetLineColor(kRed);
    l2->SetLineWidth(2);
    l2->Draw("same");
	*/

/*
  TPolyLine3D *l1 = (TPolyLine3D*) file1->Get("Line1_110");   
   l1->Draw(" same ");	
   TPolyLine3D *l2 = (TPolyLine3D*) file1->Get("Line2_110");
   l2->SetLineColor(kRed);
   l2->Draw(" same ");
  
*/
	
	/*
   TPolyLine3D *l4 = (TPolyLine3D*) file1->Get("Line4_1521");
    l4->SetLineColor(kBlue);
   l4->Draw(" same ");
 	*/
}
