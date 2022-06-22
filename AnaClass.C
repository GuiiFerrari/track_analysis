#define AnaClass_cxx
#include "AnaClass.h"


TChain* MakeChain(const TString& str2);

int main(int argc, char** argv) {

  TString str(argv[1]);

  TChain *chain = MakeChain( str);

  AnaClass t(chain);

  //t.ReadConfigFile(str);
  t.Loop(str);

  return 0;
}


TChain* MakeChain(const TString& str2) {
  auto *chain = new TChain("tree");



   //TString PathToFiles = "/home/juan/proyectos/ND_2019/pAT-TPC_17F_Sep2019/hits/";
   //TString PathToFiles = "/home/juan/proyectos/ND_2019/pAT-TPC_14O_Sep2019/hits/digi_14O/";
   TString PathToFiles = "";

  // 8B Runs
  chain->Add(PathToFiles + str2);
   //cout<<str2<<endl;



  return chain;
}




void AnaClass::ReadConfigFile(const TString& filename){

    //----get run number
    std::string temp;
    int number=0;

	for (unsigned int i=0; i < filename.Sizeof(); i++)
    		{
        //iterate the string to find the first "number" character
        //if found create another loop to extract it
        //and then break the current one
        //thus extracting the FIRST encountered numeric block
        	if (isdigit(filename[i]))
        	{
            		for (unsigned int a=i; a<filename.Sizeof(); a++)
            		{
		  		if (!isdigit(filename[a])) break;
                  		temp += filename[a];
            		}
          		break;
        		}
		}

   //cout<<temp<<endl;
   std::stringstream ss(temp);
   ss >>number;
   //cout<<"Run:  "<<number<<endl;


   //cout<<"Reading Configuration File"<<endl;
   TEnv* CFile;
   CFile = new TEnv("config/parTex.txt");


   gaspressure = CFile->GetValue("GasPressure",-1000.0);
   driftvelocity = CFile->GetValue("DriftVelocity",-1000.0);
   //nominal value ~11.46 cm (field cage)
   poszoffset = CFile->GetValue("PosZoffset",-1000.0);
   eionize = CFile->GetValue("EIonize",-1000.0);
   diffL = CFile->GetValue("DiffL",-1000.0);
   diffT = CFile->GetValue("DiffT",-1000.0);
   gain = CFile->GetValue("Gain",-1000.0);
   maxrange = CFile->GetValue("MaxRange",-1000.0);
   //40 ns sampling time
   timeperbin = CFile->GetValue("TimePerBin",-1000.0);
   Ecut = CFile->GetValue("E_PID",-1000.0);
   a0 = CFile->GetValue("A0",-1000.0);
   a1 = CFile->GetValue("A1",-1000.0);
   a2 = CFile->GetValue("A2",-1000.0);
   beam_ekin = CFile->GetValue("Beam_ekin",-1000.0);
   edead = CFile->GetValue("Edead",-1000.0);



   delete CFile;

}


void AnaClass::Loop(const TString& st)
{

 //----------------------create a new outfile for the respective run number
 string infile(st);
 string outfile("Histo" + infile);
 cout<<"Input file  "<<infile<<endl;
 cout<<"Output file  "<<outfile<<endl;
 //cout<<infile<<endl;
 cout<<"Running Analysis..."<<endl;

 file = new TFile(outfile.c_str(), "recreate");
  //cout<<infile<<endl;
  int ivt;
  float range_p1;
  vector<Float_t> vertexMean;

  TTree *anatree = new TTree("anatree","new TTree");

	anatree->Branch("ivt",&ivt);
	anatree->Branch("range_p1",&range_p1);
  anatree->Branch("vertex",&vertexMean);

  TH1F* htac1 = new TH1F("tacT","tacT",512,0,0);
  TH1F* htac2 = new TH1F("tacA","tacA",2014,0,0);
  TH1F* hprot1 = new TH1F("TriggT","TriggT",512,0,0);
  TH1F* hprot2 = new TH1F("TriggA","TriggA",512,0,0);
  TH1F* hIC1 = new TH1F("ICT","ICT",512,0,0);
  TH1F* hIC2 = new TH1F("ICA","ICA",512,0,0);
  TH1F* hmesh1 = new TH1F("meshT","meshT",512,0,0);
  TH1F* hmesh2 = new TH1F("meshA","meshA",512,0,0);
  TH1F* halpha1 = new TH1F("alphaT","alphaT",512,0,0);
  TH1F* halpha2 = new TH1F("alphaA","alphaA",512,0,0);
  TH1F* hctot = new TH1F("ctot","ctot",512,0,0);
  TH1F* hctot2 = new TH1F("ctot2","ctot2",512,0,0);
  TH2F* hpid = new TH2F("pid","pid",512,0,0,512,0,0);
  TH2F* hpid2 = new TH2F("pid2","pid2",512,0,0,512,0,0);
  TH2F* hpid3 = new TH2F("pid3","pid3",512,0,0,512,0,0);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();



   //Float_t Ecut = 1000.0;
   Float_t min_dist = 3.0; //minimum distance between point and fit line

   //Ransac = new MSimpleRansac(134,128,0);

	int contador = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //for (Long64_t jentry=0; jentry<1000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      ivt = jentry;
	

	if(TAC_A->size()>0 && TAC_A->at(0)>500 || Versor_x->size()>0) continue;
	//if(TAC_A->size()>0 && TAC_A->at(0)>500 ) continue;
      //cout<<jentry<<"  "<<TAC_A->size()<<"  "<<TAC_A->at(0)<<"  "<<vpx->size()<<"  "<<Versor_x->size()<<endl;
	cout<<jentry<<"  "<<TAC_A->size()<<"  "<<vpx->size()<<"  "<<Versor_x->size()<<endl;

	contador++;

       

	       g2D  = new TGraph2D();



	//Improving the resolution of X position (beam) by averangng the charge deposited on pixels
		int flag2 = 0;

		/*vX.clear();
		vY.clear();
		vZ.clear();
		vQ.clear();
		*/
    //cout<<vpe->size()<<"  ********************"<<endl;
	for(int j=0; j<vpx->size(); j++ ){
			int flag1 = 0;

      		if(vpe->at(j) > 110){

        /*double xcoor = vpx->at(j)*2.0; //mm
        double ycoor = 128.0 -vpy->at(j)*2.0  ; //mm
        double zoffset = 440.0; //mm
        double dvelocity = 38.0; //mm/us
        double timebin = 20.0*1e-3; //us
        double zcoor = zoffset - vpz->at(j)*dvelocity*timebin; //mm
	*/

		  /*vX.push_back(xcoor);
		  vY.push_back(ycoor);
		  vZ.push_back(vpz->at(j));
		  vQ.push_back(vpe->at(j));
			*/
      		  g2D->SetName(Form("H3D_%d",jentry));
		  g2D->SetPoint(flag2, vpx->at(j), vpt->at(j), vpy->at(j));
      		//g2D->SetPoint(flag2, xcoor, ycoor, vpz->at(j));
		  //g2D->SetPointError(flag2, 0.1, 0.1, 0.1);
		  flag2++;

      //cout<<j<<"  "<<vX[j]<<"  "<<vY[j]<<"  "<<vZ[j]<<"  "<<vpe->at(j)<<endl;
      			}
		}// for hits

		
		g2D->SetPoint(flag2, -125, 0, -125);		
		g2D->SetPoint(flag2+1, +125, 0, +125);		
		g2D->SetPoint(flag2+2, +125, 500, +125);
		
		
		
	int n = 500;
        double dt = 100.;
        //if(Nclusters<2)cout<<jentry<<"  "<<Nclusters<<"  "<<AllTracks.size()<<endl;
        lines[0]  = new TPolyLine3D(2*n);
    	lines[1]  = new TPolyLine3D(2*n);
    	lines[2]  = new TPolyLine3D(2*n);
    	lines[3]  = new TPolyLine3D(2*n);

      TVectorD integer(1);


      int flagcount = 0;
      //for(int k =0 ; k<AllTracks.size(); k++){
      int ww =0;
      for (int k = 0; k < Versor_x->size(); ++k){
        if(k>3) break;
	TVector3 newversor(Versor_x->at(k),Versor_z->at(k),Versor_y->at(k));
	TVector3 newpoint(Point_x->at(k),Point_z->at(k),Point_y->at(k));
         
        for (int i = -n; i <n;++i) {
          double t = -10 + dt*i/n;
	   
          TVector3 Linepoints =  newpoint + newversor*t;
          lines[k]->SetPoint(i,Linepoints.X(),Linepoints.Y(),Linepoints.Z());
        }
        flagcount++;
        lines[k]->Write(Form("Line_%d_%d",jentry,k));
      }

 

      anatree->Fill();
      g2D->Write();
      delete g2D;

      delete lines[0];
      delete lines[1];
      delete lines[2];
      delete lines[3];
      integer[0] = flagcount;
      //integer.Write(Form("Nlines_%d",jentry));
   }// for all events


   anatree->Write();


   htac1->Write();
   htac2->Write();
   hIC1->Write();
   hIC2->Write();
   hmesh1->Write();
   hmesh2->Write();
   hprot1->Write();
   hprot2->Write();
   halpha1->Write();
   halpha2->Write();
   hctot->Write();
   hctot2->Write();
   hpid->Write();
   hpid2->Write();
   hpid3->Write();
   file->Close();

	cout<<contador<<endl;





}
