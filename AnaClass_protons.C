#define AnaClass_protons_cxx
#include "AnaClass_protons.h"

TChain *MakeChain(const TString &str2);

int main(int argc, char **argv)
{

	TString str(argv[1]);

	TChain *chain = MakeChain(str);

	AnaClass_protons t(chain);

	t.ReadConfigFile(str);
	t.Loop(str);

	return 0;
}

TChain *MakeChain(const TString &str2)
{
	auto *chain = new TChain("tree");

	// TString PathToFiles = "/home/juan/proyectos/ND_2019/pAT-TPC_17F_Sep2019/hits/";
	//    TString PathToFiles = "/home/juan/proyectos/ND_2019/guilherme_results/data/digi_17F_plus/";
	// TString PathToFiles = "/home/juan/proyectos/ND_2019/guilherme_results/data/digi_17F_plus/";
	TString PathToFiles = "/home/ferrari/work/data/";

	// 8B Runs
	chain->Add(PathToFiles + str2);
	// cout<<str2<<endl;

	return chain;
}

void AnaClass_protons::SetERtable()
{												  // fit of the GEANT4 E vs R obtained from the simulation with the function model given by LISE++
	ifstream fER("range_17f_in_4He_350torr.txt"); // from LISE++,

	Double_t l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0;
	Int_t model = 1;
	vector<vector<Double_t>> Energy_Range;

	// LiSE++
	for (string line; getline(fER, line);)
	{
		stringstream parse_die(line);
		vector<Double_t> iRE;
		parse_die >> l1 >> l2 >> l3;
		iRE.push_back(l1); // E in MeV
		iRE.push_back(l2); // range in mm, model 1
		iRE.push_back(l3); // range in mm, model 2
		Energy_Range.push_back(iRE);
	}

	fER.close();
	Int_t v_size = Energy_Range.size();
	Double_t X[v_size];
	Double_t Y[v_size];
	for (Int_t i = 0; i < v_size; i++)
	{
		X[i] = Energy_Range.at(i).at(0);
		Y[i] = Energy_Range.at(i).at(model); //*0.738 for LISE++ Eloss to match with GEANT4
											 // cout<<X[i]<<" "<<Y[i]<<endl;
	}
	tagraph = new TGraph(v_size, Y, X);
}

void AnaClass_protons::ReadConfigFile(const TString &filename)
{

	//----get run number
	std::string temp;
	int number = 0;

	for (unsigned int i = 20; i < filename.Sizeof(); i++)
	{
		// iterate the string to find the first "number" character
		// if found create another loop to extract it
		// and then break the current one
		// thus extracting the FIRST encountered numeric block
		if (isdigit(filename[i]))
		{
			for (unsigned int a = i; a < filename.Sizeof(); a++)
			{
				if (!isdigit(filename[a]))
					break;
				temp += filename[a];
			}
			break;
		}
	}

	// cout<<temp<<endl;
	std::stringstream ss(temp);
	ss >> number;
	// cout<<"Run:  "<<number<<endl;

	// cout<<"Reading Configuration File"<<endl;
	TEnv *CFile;
	// CFile = new TEnv("config/parTex.txt");
	CFile = new TEnv(TString::Format("config/parTex_%d.par", number));

	gaspressure = CFile->GetValue("GasPressure", -1000.0);
	driftvelocity = CFile->GetValue("DriftVelocity", -1000.0);
	// nominal value ~11.46 cm (field cage)
	poszoffset = CFile->GetValue("PosZoffset", -1000.0);
	eionize = CFile->GetValue("EIonize", -1000.0);
	diffL = CFile->GetValue("DiffL", -1000.0);
	diffT = CFile->GetValue("DiffT", -1000.0);
	gain = CFile->GetValue("Gain", -1000.0);
	maxrange = CFile->GetValue("MaxRange", -1000.0);
	// 40 ns sampling time
	timeperbin = CFile->GetValue("TimePerBin", -1000.0);
	Ecut = CFile->GetValue("E_PID", -1000.0);
	a0 = CFile->GetValue("A0", -1000.0);
	a1 = CFile->GetValue("A1", -1000.0);
	a2 = CFile->GetValue("A2", -1000.0);
	beam_ekin = CFile->GetValue("Beam_ekin", -1000.0);
	edead = CFile->GetValue("Edead", -1000.0);
	tb0 = CFile->GetValue("TB0", -1000.0);

	delete CFile;
}

double AnaClass_protons::PointLineDist(int i, TVector3 Vs, TVector3 Ps)
{
	// distance point to line
	TVector3 newPoint = {vpx->at(i), vpy->at(i), (tb0 - vpt->at(i)) * driftvelocity};
	TVector3 vec = Ps - newPoint;
	TVector3 nD = Vs.Cross(vec);
	double dist = nD.Mag() / Vs.Mag();

	return dist;
}

void AnaClass_protons::Loop(const TString &st)
{

	//----------------------create a new outfile for the respective run number
	string infile(st);
	string outfile("Histo" + infile);
	// cout<<"Input file  "<<infile<<endl;
	// cout<<"Output file  "<<outfile<<endl;
	// cout<<infile<<endl;
	// cout<<"Running Analysis..."<<endl;

	TFile *fcut = new TFile("cutpid2.root");
	TCutG *CUTCharge = (TCutG *)fcut->Get("cutp2");
	// TCutG *CUTCharge = (TCutG *) fcut->Get("cuta2");
	// TCutG *CUTCharge = (TCutG *) fcut->Get("cutu");

	TFile *fcut2 = new TFile("cutbeam.root");
	TCutG *CUTBeam = (TCutG *)fcut2->Get("cutbeam");

	TFile *fcut3 = new TFile("cutbeamCh3.root");
	TCutG *CUTBeamCh = (TCutG *)fcut3->Get("cutbeamCh3");

	TFile *fcut4 = new TFile("cutbeamCh4.root");
	TCutG *CUTBeamCh4 = (TCutG *)fcut4->Get("cutbeamCh4");

	file = new TFile(outfile.c_str(), "recreate");

	file->mkdir("Th_zprof");
	// file->mkdir("Th_zprof_coin");
	file->mkdir("Zvert_thprof");

	// cout<<infile<<endl;
	int ivt;
	float range_p1;
	vector<Float_t> vertexMean;

	TTree *anatree = new TTree("anatree", "new TTree");

	anatree->Branch("ivt", &ivt);
	anatree->Branch("range_p1", &range_p1);
	anatree->Branch("vertex", &vertexMean);

	TH1F *htac1 = new TH1F("tacT", "tacT", 512, 0, 0);
	TH1F *htac2 = new TH1F("tacA", "tacA", 2014, 0, 0);
	TH1F *hprot1 = new TH1F("TriggT", "TriggT", 512, 0, 0);
	TH1F *hprot2 = new TH1F("TriggA", "TriggA", 512, 0, 0);
	TH1F *hIC1 = new TH1F("ICT", "ICT", 512, 0, 0);
	TH1F *hIC2 = new TH1F("ICA", "ICA", 512, 0, 0);
	TH1F *hmesh1 = new TH1F("meshT", "meshT", 512, 0, 0);
	TH1F *hmesh2 = new TH1F("meshA", "meshA", 512, 0, 0);
	TH1F *halpha1 = new TH1F("alphaT", "alphaT", 512, 0, 0);
	TH1F *halpha2 = new TH1F("alphaA", "alphaA", 512, 0, 0);
	TH1F *hctot = new TH1F("ctot", "ctot", 512, 0, 0);
	TH1F *hctot2 = new TH1F("ctot2", "ctot2", 512, 0, 0);
	TH2F *hpid = new TH2F("pid", "pid", 512, 0, 300, 512, 0, 16e3);
	TH2F *hpid2 = new TH2F("pid2", "pid2", 512, 0, 300, 512, 0, 5e3);
	TH2F *hpid3 = new TH2F("pid3", "pid3", 512, 0, 5e3, 512, 0, 16e3);
	TH1F *hunrea = new TH1F("unrea", "unrea beam", 300, -200, 700);
	TH1F *hverZ = new TH1F("verZ", "verZ", 300, -200, 800);
	TH1F *hlastZ = new TH1F("lastZ", "lastZ", 300, -200, 800);
	TH1F *hdE = new TH1F("hdE", "hdE", 300, 0, 0);
	TH1F *hTheta = new TH1F("hTheta", "hTheta", 180, 0, 180);
	TH1F *hLength = new TH1F("hLength", "hLength", 180, 0, 0);
	TH2F *hTvsL = new TH2F("Comprimento (mm) vs Angulo (graus)",
						   "Comprimento (mm) vs Angulo (graus)", 180, 0, 180, 512, 0, 500);
	hTvsL->GetXaxis()->SetTitle("Angulo (graus)");
	hTvsL->GetYaxis()->SetTitle("Comprimento (mm)");
	TH2F *hTvsLCoin = new TH2F("Comprimento (mm) vs Angulo (graus) (coinc)",
							   "Comprimento (mm) vs Angulo (graus)", 180, 0, 180, 512, 0, 500);
	hTvsLCoin->GetXaxis()->SetTitle("Angulo (graus)");
	hTvsLCoin->GetYaxis()->SetTitle("Comprimento (mm)");
	TH2F *hTvsCh = new TH2F("Carga vs Angulo (graus)",
							"Carga vs Angulo (graus)", 180, 0, 180, 512, 0, 4e3);
	hTvsCh->GetXaxis()->SetTitle("Angulo (graus)");
	hTvsCh->GetYaxis()->SetTitle("Carga");
	TH2F *hLvsCh_hit = new TH2F("LvsCh_hit", "LvsCh_hit", 300, 0, 500, 512, 0, 2000);
	TH2F *hLvsCh = new TH2F("Carga vs Comprimento (mm)",
							"Carga vs Comprimento (mm)", 300, 0, 500, 512, 0, 2.5e5);
	hLvsCh->GetXaxis()->SetTitle("Comprimento (mm)");
	hLvsCh->GetYaxis()->SetTitle("Carga");
	TH2F *hChvsAlpha = new TH2F("ChvsAlpha", "ChvsAlpha", 500, 0, 16e4, 512, 0, 6000);
	TH1F *hNliers_eve = new TH1F("Nliers_eve", "Nliers_eve", 512, 0, 0);
	TH1F *hChproj = new TH1F("Chproj", "Chproj", 500, 0, 6e4);
	TH2F *hTh_Zvert = new TH2F("Th_Zvert", "Th_Zvert", 300, -50, 600, 180, 0, 180);
	TH1F *hEprotons = new TH1F("En_protons", "En_protons", 50, 0, 500);
	TH1F *hNTracks = new TH1F("NTracks", "NTracks", 20, 0, 20);

	for (int i = 0; i < 50; i++)
	{
		int w = i * 10;
		thetaHisto[i] = new TH1F(Form("Theta_hist_%d", w), "", 45, 0, 180);
		EBeamHisto[i] = new TH1F(Form("EBeam_hist_%d", w), "", 200, -10, 50);
		thetaHisto_coin[i] = new TH1F(Form("Theta_hist_C_%d", w), "", 45, 0, 180);
	}
	for (int i = 0; i < 20; i++)
	{
		int w = i * 5;
		ZverHisto[i] = new TH1F(Form("Zvert_hist_%d", w), "", 500, 0, 500);
	}

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntriesFast();

	SetERtable();

	// Float_t Ecut = 1000.0;
	Float_t min_dist = 3.0;	  // minimum distance between point and fit line
	double luminosity = 2.49; // 1/mb

	// Ransac = new MSimpleRansac(134,128,0);

	int contador = 0;
	double thetaMin = 0;
	double step_theta = 4.0; // in degrees
	double step_thetaInt = 5.0;
	double rangeMin = 0;
	double step_range = 10.0; // 20 cm
	double mindist_vert = 20.0;

	// driftvelocity = 1.0;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		// for (Long64_t jentry=0; jentry<1000;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		ivt = jentry;

		if (TAC_A->size() > 0 && TAC_A->at(0) > 500 || Versor_x->size() == 0)
			continue; // solamente eventos de reaccion
		// if(IC_T->at(0)>90 || IC_T->at(0)<70 ) continue;
		// if(TAC_A->size()>0 && TAC_A->at(0)<500) continue;
		// if(TAC_A->size()>0 && TAC_A->at(0)>500) continue;
		// cout<<jentry<<"  "<<TAC_A->size()<<"  "<<TAC_A->at(0)<<"  "<<vpx->size()<<"  "<<Versor_x->size()<<endl;
		// cout<<jentry<<"  "<<TAC_A->size()<<"  "<<vpx->size()<<"  "<<Versor_x->size()<<endl;

		contador++;
		double dEcharge = 0;

		for (int j = 0; j < vpx->size(); j++)
		{

			double coorZ = (tb0 - vpt->at(j)) * driftvelocity;
			// cout<<coorZ<<"  "<<driftvelocity<<endl;
			if (sqrt(vpx->at(j) * vpx->at(j) + vpy->at(j) * vpy->at(j)) < 30)
			{
				hunrea->Fill(coorZ, vpe->at(j));
				// hunrea->Fill(vpt->at(j),vpe->at(j));
				if (coorZ > 40 && coorZ < 100)
					dEcharge += vpe->at(j);
			}
		}

		hdE->Fill(dEcharge);

		double gainA = 1;

		for (int j = 0; j < TAC_T->size(); j++)
		{
			if (Trigg_T->size() != 1)
				continue;
			if (Trigg_T->at(0) < 400)
				continue;
			if (TAC_T->at(0) < 400)
				continue;
			if (TAC_T->size() != 1)
				continue;
			if (Mesh_T->size() != 1)
				continue;
			// if(Trigg_T->size()!=1) continue;
			// if(alpha_A->at(0)==0) gainA = 1;
			// else gainA = alpha_A->at(0)*2/1000.0;
			hmesh1->Fill(Mesh_T->at(j));
			hmesh2->Fill(Mesh_A->at(j) * gainA);
			htac1->Fill(TAC_T->at(j));
			htac2->Fill(TAC_A->at(j));
			hpid->Fill(TAC_A->at(j), dEcharge);
			hpid2->Fill(TAC_A->at(j), Mesh_A->at(j) * gainA);
			hpid3->Fill(Mesh_A->at(j), dEcharge);
			hprot1->Fill(Trigg_T->at(j));
			hprot2->Fill(Trigg_A->at(j));
			// halpha1->Fill(alpha_T->at(0));
			// halpha2->Fill(alpha_A->at(0));
		}

		for (int j = 0; j < IC_T->size(); j++)
		{
			hIC1->Fill(IC_T->at(j));
			hIC2->Fill(IC_A->at(j));
			// hdE->Fill(IC_A->at(j));
			// hpid->Fill(TAC_A->at(j),dEcharge);
		}

		// cout<<TAC_T->size()<<"  "<<Trigg_T->size()<<"  "<<Mesh_T->size()<<"  "<<alpha_A->at(0)<<endl;

		TVector3 BeamDir(0., 0., 1.0);
		double theta = 0;
		double totalCharge = 0;
		std::vector<bool> flagTrack;
		flagTrack.clear();
		for (size_t i = 0; i < Vertex_z->size(); i++)
		{
			flagTrack.push_back(false);
		}

		for (int j = 0; j < Vertex_z->size(); j++)
		{
			double vertexpointZ1 = (tb0 - Vertex_z->at(j)) * driftvelocity;
			TVector3 theVertex1(Vertex_x->at(j), Vertex_y->at(j), vertexpointZ1);
			// double nli_charge1 = Num_inliers->at(j)*Total_charge->at(j);
			double nli_charge1 = Total_charge->at(j);
			// if(vertexpointZ1<40 || vertexpointZ1>500) continue;
			if (Num_inliers->at(j) < 30)
				continue;
			if (Vertex_z->size() == 1 && (vertexpointZ1 > 0 && vertexpointZ1 < 500) && Is_primary->at(j) == 0)
				flagTrack[j] = true;
			// if(Vertex_z->size()==1  && Is_primary->at(j)==0) flagTrack[j]=true;

			for (int i = j + 1; i < Vertex_z->size(); i++)
			{
				double vertexpointZ2 = (tb0 - Vertex_z->at(i)) * driftvelocity;
				TVector3 theVertex2(Vertex_x->at(i), Vertex_y->at(i), vertexpointZ2);
				// double nli_charge2 = Num_inliers->at(i)*Total_charge->at(i);
				double nli_charge2 = Total_charge->at(i);
				// if(vertexpointZ2<40 || vertexpointZ2>500) continue;
				// if(Is_primary->at(i)==1) continue;
				TVector3 diffVert = theVertex1 - theVertex2;
				if (diffVert.Mag() > mindist_vert && nli_charge1 >= nli_charge2)
				{
					flagTrack[j] = true;
					flagTrack[i] = false;
				}
				if (diffVert.Mag() > mindist_vert && nli_charge1 < nli_charge2)
				{
					flagTrack[j] = false;
					flagTrack[i] = true;
				}
				if (diffVert.Mag() < mindist_vert)
				{
					flagTrack[j] = true;
					flagTrack[i] = true;
				}
				// if(vertexpointZ2<30 /*|| vertexpointZ2>500*/ || Is_primary->at(i)==1) flagTrack[i]=false;
				if (Is_primary->at(i) == 1)
					flagTrack[i] = false;

				// cout<<jentry<<"  "<<j<<"  "<<i<<"  "<<diffVert.Mag()<<endl;
			}
			if (vertexpointZ1 < 0 || vertexpointZ1 > 500 || Is_primary->at(j) == 1)
				flagTrack[j] = false;
			// if( Is_primary->at(j)==1) flagTrack[j]=false;
		}

		int tracks_filter1 = 0;
		for (int j = 0; j < Vertex_z->size(); j++)
		{
			if (flagTrack[j] == true)
				tracks_filter1++;
		}

		if (tracks_filter1 > 1)
		{
			// confirmar tracks con el mismo vertice
			for (int j = 0; j < Vertex_z->size(); j++)
			{
				double vertexpointZ1 = (tb0 - Vertex_z->at(j)) * driftvelocity;
				TVector3 theVertex1(Vertex_x->at(j), Vertex_y->at(j), vertexpointZ1);
				TVector3 theVersor1(Versor_x->at(j), Versor_y->at(j), Versor_z->at(j));
				double nli_charge1 = Total_charge->at(j);
				// double nli_charge1 = Num_inliers->at(j);
				if (flagTrack[j] == false)
					continue;

				// if(nli_charge1<30){ flagTrack[j]=false; continue;}

				for (int i = j + 1; i < Vertex_z->size(); i++)
				{
					double vertexpointZ2 = (tb0 - Vertex_z->at(i)) * driftvelocity;
					TVector3 theVertex2(Vertex_x->at(i), Vertex_y->at(i), vertexpointZ2);
					TVector3 theVersor2(Versor_x->at(i), Versor_y->at(i), Versor_z->at(i));
					double nli_charge2 = Total_charge->at(i);
					// double nli_charge2 = Num_inliers->at(i);
					double angulo12 = theVersor1.Angle(theVersor2) * 180. / 3.1415;
					if (flagTrack[i] == false)
						continue;
					// if(nli_charge2<30){ flagTrack[i]=false; continue;}

					TVector3 diffVert = theVertex1 - theVertex2;
					if (diffVert.Mag() > mindist_vert && nli_charge1 >= nli_charge2)
					{
						flagTrack[j] = true;
						flagTrack[i] = false;
					}
					if (diffVert.Mag() > mindist_vert && nli_charge1 < nli_charge2)
					{
						flagTrack[j] = false;
						flagTrack[i] = true;
					}
					if (diffVert.Mag() < mindist_vert && (angulo12 < 10 || angulo12 > 170))
					{
						if (nli_charge1 >= nli_charge2)
						{
							flagTrack[j] = true;
							flagTrack[i] = false;
						}
						if (nli_charge1 < nli_charge2)
						{
							flagTrack[j] = false;
							flagTrack[i] = true;
						}
					}
				}
			}
		}

		int tracks_filter2 = 0;
		for (int j = 0; j < Vertex_z->size(); j++)
		{
			if (flagTrack[j] == true)
				tracks_filter2++;
		}

		double gainAlfas = 1;
		double longi = 0;
		double carga = 0;
		int nliers = 0;

		// if(tracks_filter2>2) cout<<jentry<<"  "<<Vertex_z->size()<<"  "<<tracks_filter1<<"  "<<tracks_filter2<<endl;

		bool proton_flag = false;

		for (int j = 0; j < Vertex_z->size(); j++)
		{
			if (flagTrack[j] == false)
				continue;
			if (Trigg_T->size() != 1)
				continue;
			if (Num_inliers->at(j) < 40)
				continue;
			// if(Trigg_T->at(0)<400) continue;
			// if(alpha_A->at(0)>500) continue; //67% de la estadistica esta en alpha>500

			double vertexpointZ = (tb0 - Vertex_z->at(j)) * driftvelocity;
			if (vertexpointZ > 500 || vertexpointZ < 0)
				continue;

			double lastpointZ = (tb0 - Max_Pos_z->at(j)) * driftvelocity;
			double versorZ = -1 * (Versor_z->at(j)) * driftvelocity;
			double linepointZ = (tb0 - Point_z->at(j)) * driftvelocity;

			TVector3 theVertex(Vertex_x->at(j), Vertex_y->at(j), vertexpointZ);
			TVector3 lastP(Max_Pos_x->at(j), Max_Pos_y->at(j), lastpointZ);
			TVector3 theVersor(TMath::Sign(1, lastP.X()) * fabs(Versor_x->at(j)), TMath::Sign(1, lastP.Y()) * fabs(Versor_y->at(j)), TMath::Sign(1, (lastpointZ - vertexpointZ)) * fabs(versorZ));
			TVector3 thePoint(Point_x->at(j), Point_y->at(j), linepointZ);
			double angle = theVersor.Angle(BeamDir) * 180. / 3.1415;

			double plast_radius = sqrt(pow(lastP.X(), 2.0) + pow(lastP.Y(), 2.0));
			// if(lastP.Z()>500 || lastP.Z()<0 || plast_radius>125) continue;
			if (lastP.Z() > 500 || lastP.Z() < 0)
				continue;

			// if(alpha_T->size()>0) halpha1->Fill(alpha_T->at(0));
			// if(alpha_T->size()>0 && alpha_A->at(0)>500){ halpha2->Fill(alpha_A->at(0));  gainAlfas = alpha_A->at(0)*2/1000.0;}

			// for(int w=0; w<Trigg_T->size(); w++ ){
			// if(Trigg_T->size()==1) continue;
			TVector3 TraLength = lastP - theVertex;
			// Função abaixo verifica se o ponto (comprimento, carga * ganho)
			// está dentro da região de corte
			if (!CUTCharge->IsInside(TraLength.Mag(), gain * Total_charge->at(j)))
				continue;
			if (CUTBeamCh4->IsInside(TraLength.Mag(), gain * Total_charge->at(j) / Num_inliers->at(j)))
				continue; // quitar el beam scattered
			if (CUTBeam->IsInside(angle, TraLength.Mag()))
				continue; // quitar el beam scattered
			//////// if(CUTBeamCh->IsInside(angle,Total_charge->at(j)*gainAlfas/Num_inliers->at(j))) continue;	//quitar el beam scattered
			if (vertexpointZ > 30 && vertexpointZ < 70 && angle > 10 && angle < 100)
				continue; // quitar un alfa elastic que sale del volumen

			//////// if(angle<20) continue; //para seleccionar protones y quitar beam
			// if (angle < 5 || angle > 175)
			// 	continue; // para quitar fits paralelos al beam

			proton_flag = true;

			// Range to Ekin
			double eLoss_reco = tagraph->Eval(vertexpointZ);
			// cout<<vertexpointZ<<"  "<<eLoss_reco<<endl;
			double Ebeam = 34.76 - eLoss_reco;

			hTvsL->Fill(angle, TraLength.Mag());
			hTvsCh->Fill(angle, Total_charge->at(j) * gainAlfas / Num_inliers->at(j));
			hLength->Fill(TraLength.Mag());
			hverZ->Fill(vertexpointZ);
			hlastZ->Fill(lastpointZ);
			hLvsCh_hit->Fill(TraLength.Mag(), gain * Total_charge->at(j) / Num_inliers->at(j));
			hLvsCh->Fill(TraLength.Mag(), gain * Total_charge->at(j));
			hChvsAlpha->Fill(Total_charge->at(j), alpha_A->at(0) + 1000);
			if (TraLength.Mag() > 80 && TraLength.Mag() < 85)
				hChproj->Fill(gain * Total_charge->at(j));
			longi = TraLength.Mag();
			carga = gain * Total_charge->at(j) / Num_inliers->at(j);
			theta = angle;
			nliers = Num_inliers->at(j);
			hNliers_eve->Fill(nliers);
			hTh_Zvert->Fill(vertexpointZ, angle);

			halpha1->Fill(alpha_T->at(0));
			halpha2->Fill(alpha_A->at(0));

			int ThetaBin = TMath::Floor((angle - thetaMin) / step_theta);
			double solid_ang = 2.0 * 3.14159 * (cos(step_theta * (ThetaBin)*3.141592 / 180.0) - cos(step_theta * (ThetaBin + 1) * 3.141592 / 180.0));
			int RangeBin = TMath::Floor((vertexpointZ - rangeMin) / step_range);

			// protones
			// if(longi>100 && carga<500 && nliers>30) hTheta->Fill(angle,1.0/solid_ang);
			// if(longi>100 && carga<500 && nliers>30) thetaHisto[RangeBin]->Fill(angle,1.0/solid_ang);
			// alfas
			// if(longi<120 &&  nliers>30) hTheta->Fill(angle,1.0/solid_ang);
			// if(longi<120  && nliers>30) thetaHisto[RangeBin]->Fill(angle,1.0/solid_ang);
			hTheta->Fill(angle, 1.0 / (solid_ang * luminosity));
			// thetaHisto[RangeBin]->Fill(angle,1.0/(solid_ang*luminosity));
			thetaHisto[RangeBin]->Fill(angle);
			EBeamHisto[RangeBin]->Fill(Ebeam);

			int ThetaBinInt = TMath::Floor((angle - thetaMin) / step_thetaInt);
			if (angle < 100)
				ZverHisto[ThetaBinInt]->Fill(vertexpointZ);

			//}
		}

		int tracks_filter3 = 0;
		int coincidences = 0;

		for (int j = 0; j < Vertex_z->size(); j++)
		{
			if (proton_flag == false)
				continue; // there is a proton track in the event
			if (tracks_filter2 < 2)
				continue; // more than 2 good tracks in the event
			if (flagTrack[j] == false)
				continue; // is a good track?
			if (Num_inliers->at(j) < 30)
				continue;

			double vertexpointZ = (tb0 - Vertex_z->at(j)) * driftvelocity;
			if (vertexpointZ > 500 || vertexpointZ < 0)
				continue;

			double lastpointZ = (tb0 - Max_Pos_z->at(j)) * driftvelocity;
			double versorZ = -1 * (Versor_z->at(j)) * driftvelocity;
			double linepointZ = (tb0 - Point_z->at(j)) * driftvelocity;

			TVector3 theVertex(Vertex_x->at(j), Vertex_y->at(j), vertexpointZ);
			TVector3 lastP(Max_Pos_x->at(j), Max_Pos_y->at(j), lastpointZ);
			TVector3 theVersor(TMath::Sign(1, lastP.X()) * fabs(Versor_x->at(j)), TMath::Sign(1, lastP.Y()) * fabs(Versor_y->at(j)), TMath::Sign(1, (lastpointZ - vertexpointZ)) * fabs(versorZ));
			TVector3 thePoint(Point_x->at(j), Point_y->at(j), linepointZ);
			double angle = theVersor.Angle(BeamDir) * 180. / 3.1415;

			double plast_radius = sqrt(pow(lastP.X(), 2.0) + pow(lastP.Y(), 2.0));
			// if(lastP.Z()>500 || lastP.Z()<0 || plast_radius>125) continue;
			if (lastP.Z() > 500 || lastP.Z() < 0)
				continue;

			TVector3 TraLength = lastP - theVertex;

			// if(angle<20) continue; //para seleccionar protones y quitar beam
			if (angle < 5 || angle > 175)
				continue; // para quitar fits paralelos al beam

			for (int i = j + 1; i < Vertex_z->size(); i++)
			{
				double vertexpointZ2 = (tb0 - Vertex_z->at(i)) * driftvelocity;
				double lastpointZ2 = (tb0 - Max_Pos_z->at(i)) * driftvelocity;
				double versorZ2 = -1 * (Versor_z->at(i)) * driftvelocity;
				TVector3 theVertex2(Vertex_x->at(i), Vertex_y->at(i), vertexpointZ2);
				TVector3 lastP2(Max_Pos_x->at(i), Max_Pos_y->at(i), lastpointZ2);

				TVector3 theVersor2(TMath::Sign(1, lastP2.X()) * fabs(Versor_x->at(i)), TMath::Sign(1, lastP2.Y()) * fabs(Versor_y->at(i)), TMath::Sign(1, (lastpointZ2 - vertexpointZ2)) * fabs(versorZ2));

				double nli_charge2 = Total_charge->at(i);

				double angle2 = theVersor2.Angle(BeamDir) * 180. / 3.1415;
				TVector3 TraLength2 = lastP2 - theVertex2;
				// double nli_charge2 = Num_inliers->at(i);
				double angulo12 = theVersor.Angle(theVersor2) * 180. / 3.1415;

				TVector3 diffVert = theVertex - theVertex2;

				if (diffVert.Mag() < mindist_vert)
				{
					if (CUTBeam->IsInside(angle, TraLength.Mag()) || CUTBeam->IsInside(angle2, TraLength2.Mag()))
						coincidences++;
					// one of the particles is a beam
				}
			}

			// Range to Ekin
			double eLoss_reco = tagraph->Eval(vertexpointZ);
			// cout<<vertexpointZ<<"  "<<eLoss_reco<<endl;
			double Ebeam = 34.76 - eLoss_reco;

			if (coincidences == 1)
				hTvsLCoin->Fill(angle, TraLength.Mag());
			tracks_filter3 = coincidences + 1;

			int ThetaBin = TMath::Floor((angle - thetaMin) / step_theta);
			double solid_ang = 2.0 * 3.14159 * (cos(step_theta * (ThetaBin)*3.141592 / 180.0) - cos(step_theta * (ThetaBin + 1) * 3.141592 / 180.0));
			int RangeBin = TMath::Floor((vertexpointZ - rangeMin) / step_range);

			// if(!CUTBeam->IsInside(angle,TraLength.Mag()) && coincidences==3) thetaHisto_coin[RangeBin]->Fill(angle,1.0/(solid_ang*luminosity));
			if (!CUTBeam->IsInside(angle, TraLength.Mag()) && coincidences == 1)
				thetaHisto_coin[RangeBin]->Fill(angle);
		}

		if (tracks_filter3 > 0)
			hNTracks->Fill(tracks_filter3);

		/*
			TVector3 TraLength = lastP - theVertex;
			hTvsL->Fill(angle,TraLength.Mag());
			hTvsCh->Fill(angle,Total_charge->at(j)/Num_inliers->at(j));
			hTheta->Fill(angle);
			hLength->Fill(TraLength.Mag());
			hverZ->Fill(vertexpointZ);
			hlastZ->Fill(lastpointZ);
			theta = angle;
			//break;
			*/

		// cout<<jentry<<"  "<<Vertex_z->size()<<"  "<<counter<<endl;

		/*
		if(theta<20 && nliers>30 ){
			cout<<jentry<<"  "<<theta<<"  "<<Trigg_T->size()<<endl;
//-----------------------------------crear figuras 3D-------------------------------------------------
		   g2D  = new TGraph2D();
		   g2D->SetName(Form("H3D_%d",jentry));
		   g2D_vertex = new TGraph2D();
		   g2D_vertex->SetName(Form("vertex_%d",jentry));
		   g2D_lpoint  = new TGraph2D();
		   g2D_lpoint->SetName(Form("lpoint_%d",jentry));


		int flag2 = 0;

	for(int j=0; j<vpx->size(); j++ ){
			int flag1 = 0;
		double coorZ = (tb0-vpt->at(j))*driftvelocity;

			if(vpe->at(j) > 110){


		  //g2D->SetPoint(flag2, vpx->at(j), vpt->at(j), vpy->at(j));
		  g2D->SetPoint(flag2, vpx->at(j), coorZ, vpy->at(j));
		  flag2++;

				}
		}// for hits


		g2D->SetPoint(flag2, -125, 0, -125);
		g2D->SetPoint(flag2+1, +125, 0, +125);
		g2D->SetPoint(flag2+2, +125, 600, +125);



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

	TVector3 newversor(Versor_x->at(k),-1*(Versor_z->at(k))*driftvelocity,Versor_y->at(k));
	TVector3 newpoint(Point_x->at(k),(tb0-Point_z->at(k))*driftvelocity,Point_y->at(k));
	//TVector3 newversor(Versor_x->at(k),Versor_z->at(k),Versor_y->at(k));
	//TVector3 newpoint(Point_x->at(k),Point_z->at(k),Point_y->at(k));


	g2D_vertex->SetPoint(k, Vertex_x->at(k), (tb0-Vertex_z->at(k))*driftvelocity, Vertex_y->at(k));
	g2D_lpoint->SetPoint(k, Max_Pos_x->at(k), (tb0-Max_Pos_z->at(k))*driftvelocity, Max_Pos_y->at(k));
	//g2D_vertex->SetPoint(k, Vertex_x->at(k), Vertex_z->at(k), Vertex_y->at(k));
	//g2D_lpoint->SetPoint(k, Max_Pos_x->at(k), Max_Pos_z->at(k), Max_Pos_y->at(k));
	//cout<<jentry<<"  "<<k<<"  "<<Vertex_x->at(k)<<"  "<<Vertex_y->at(k)<<"  "<<Vertex_z->at(k)<<endl;
	//cout<<jentry<<"  "<<k<<"  "<<Max_Pos_x->at(k)<<"  "<<Max_Pos_y->at(k)<<"  "<<Max_Pos_z->at(k)<<endl;

		for (int i = -n; i <n;++i) {
		  double t = -10 + dt*i/n;

		  TVector3 Linepoints =  newpoint + newversor*t;
	  if(flagTrack[k]==true)
		  lines[k]->SetPoint(i,Linepoints.X(),Linepoints.Y(),Linepoints.Z());
		}
	//if(flagTrack[k]==false) continue;
		flagcount++;
		lines[k]->Write(Form("Line_%d_%d",jentry,k));
	  }

	g2D_vertex->SetPoint(flagcount, 0, 500, 0);
	g2D_lpoint->SetPoint(flagcount, 0, 500, 0);

	  anatree->Fill();
	  g2D->Write();
	  delete g2D;
	  g2D_vertex->Write();
	  delete g2D_vertex;
	  g2D_lpoint->Write();
	  delete g2D_lpoint;

	  delete lines[0];
	  delete lines[1];
	  delete lines[2];
	  delete lines[3];
	  integer[0] = flagcount;
	  integer.Write(Form("Nlines_%d",jentry));
///-----------------------------------------end of crear figuras 3D
	}
		*/
	} // for all events

	/*
	hunrea->GetXaxis()->SetRangeUser(0, 450);
	int binmax = hunrea->GetMaximumBin();
	double x = hunrea->GetXaxis()->GetBinCenter(binmax);
	hunrea->GetXaxis()->SetRangeUser(240, 310);
	TF1 *f1 = new TF1("f1","gaus",240,310);
	hunrea->Fit("f1", "Q0");

	cout<<infile<<"  "<<x<<"  "<<f1->GetParameter(1)<<"  "<<f1->GetParameter(2)<<endl;
	*/
	cout << infile << "  " << ivt << "  " << nentries << endl;

	for (int i = 0; i < 50; i++)
	{

		int w = i * 10;
		double eman = EBeamHisto[i]->GetMean();
		double sigmaE = EBeamHisto[i]->GetRMS();
		// double countsBinE = EBeamHisto[i]->GetEntries()/(2.0*sigmaE);
		double countsBinE = EBeamHisto[i]->GetEntries();
		// if(countsBinE>0)
		// hEprotons->Fill(w,countsBinE/luminosity);
		hEprotons->Fill(w, countsBinE);
	}

	anatree->Write();

	hTvsCh->Write();
	hTvsL->Write();
	hTvsLCoin->Write();
	hNTracks->Write();
	hLength->Write();
	hTheta->Write();
	hunrea->Write();
	hverZ->Write();
	hlastZ->Write();
	hdE->Write();
	hEprotons->Write();
	hLvsCh_hit->Write();
	hLvsCh->Write();
	hChvsAlpha->Write();
	hChproj->Write();
	hNliers_eve->Write();
	hTh_Zvert->Write();
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
	file->cd("Th_zprof");
	for (int i = 0; i < 50; i++)
	{
		thetaHisto[i]->Write();
		EBeamHisto[i]->Write();
		thetaHisto_coin[i]->Write();
	}
	file->cd("Zvert_thprof");
	for (int i = 0; i < 20; i++)
		ZverHisto[i]->Write();
	file->Close();

	// cout<<contador<<endl;
}
