//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 29 14:08:00 2021 by ROOT version 6.12/06
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef AnaClass2_h
#define AnaClass2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>



#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include "TF1.h"
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TH2Poly.h>
#include "TEnv.h"
#include <TRandom3.h>


#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include <iomanip>
#include <sstream>
#include <math.h>

using namespace std;

class AnaClass2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *Event_number;
   vector<double>  *Versor_x;
   vector<double>  *Versor_y;
   vector<double>  *Versor_z;
   vector<double>  *Point_x;
   vector<double>  *Point_y;
   vector<double>  *Point_z;
   vector<double>  *Vertex_x;
   vector<double>  *Vertex_y;
   vector<double>  *Vertex_z;
   vector<int>     *Is_primary;
   vector<double>  *Total_charge;
   vector<double>  *Max_Pos_x;
   vector<double>  *Max_Pos_y;
   vector<double>  *Max_Pos_z;
   vector<int>     *Num_inliers;
   vector<double>  *vpx;
   vector<double>  *vpy;
   vector<double>  *vpz;
   vector<double>  *vpt;
   vector<double>  *vpe;
   vector<double>  *vhits;
   vector<double>  *TAC_T;
   vector<double>  *TAC_A;
   vector<double>  *Trigg_T;
   vector<double>  *Trigg_A;
   vector<double>  *IC_T;
   vector<double>  *IC_A;
   vector<double>  *Mesh_T;
   vector<double>  *Mesh_A;
   vector<double>  *alpha_T;
   vector<double>  *alpha_A;

   // List of branches
   TBranch        *b_Event_number;   //!
   TBranch        *b_Versor_x;   //!
   TBranch        *b_Versor_y;   //!
   TBranch        *b_Versor_z;   //!
   TBranch        *b_Point_x;   //!
   TBranch        *b_Point_y;   //!
   TBranch        *b_Point_z;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Is_primary;   //!
   TBranch        *b_Total_charge;   //!
   TBranch        *b_Max_Pos_x;   //!
   TBranch        *b_Max_Pos_y;   //!
   TBranch        *b_Max_Pos_z;   //!
   TBranch        *b_Num_inliers;   //!
   TBranch        *b_vpx;   //!
   TBranch        *b_vpy;   //!
   TBranch        *b_vpz;   //!
   TBranch        *b_vpt;   //!
   TBranch        *b_vpe;   //!
   TBranch        *b_vhits;   //!
   TBranch        *b_TAC_T;   //!
   TBranch        *b_TAC_A;   //!
   TBranch        *b_Trigg_T;   //!
   TBranch        *b_Trigg_A;   //!
   TBranch        *b_IC_T;   //!
   TBranch        *b_IC_A;   //!
   TBranch        *b_Mesh_T;   //!
   TBranch        *b_Mesh_A;   //!
   TBranch        *b_alpha_T;   //!
   TBranch        *b_alpha_A;   //!

   AnaClass2(TTree *tree=0);
   virtual ~AnaClass2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const TString& st);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void ReadConfigFile(const TString& filename);
   double PointLineDist(int i,TVector3 Vs, TVector3 Ps);
   void SetERtable();
   void getkin(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double Ex4);
   double sq(double val );
   double omega(double x, double y, double z);
   void TwobodyKin(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject);
   double getThetaLab(double m1, double m2, double m3, double m4, double K_proj, double thetacm, double Ex);
   double getEkin2(double m1, double m2, double m3, double m4, double K_proj, double thetacm, double Ex);

   TFile *file;

   Float_t gaspressure;
   Float_t driftvelocity;
   Float_t poszoffset;
   Float_t eionize;
   Float_t diffL;
   Float_t diffT;
   Float_t gain;
   Float_t maxrange;
   Float_t timeperbin;
   Float_t Ecut;
   Float_t a0;
   Float_t a1;
   Float_t a2;
   Float_t beam_ekin;
   Float_t edead;
   Float_t tb0;
   TGraph *tagraph;
   TGraph *tagraph_4He;
   double ThetaCM;
   double EkinScatt;
   double Ex_reco;

   vector<double> vX;
   vector<double> vY;
   vector<double> vZ;
   vector<double> vZ_smooth;
   vector<double> vQ;
   vector<double> vQ_buffer;
   vector<double> vX_buffer;
   vector<double> vY_buffer;
   vector<double> vZ_buffer;

   std::vector< std::vector<float> *> pts;
   std::vector<float>* p;
   std::vector<std::vector<float> *> * mModels;
   std::vector<unsigned int> Lables;
   std::vector<unsigned int> LableCount;

   struct Track3d{
     //vector<int> inIndx;
     TVector3 Vl;
     TVector3 Pl;
     TVector3 Vertex;	
     double length;
     double charge;
     double Theta;
   };

   typedef std::vector<Track3d> Track3Dvector;
   Track3Dvector AllTracks;


   vector<bool> vBeam;
   vector<double> vQTtrack;

   TVector3 Vdirec;
   TPolyLine3D * lines[4];

   vector<double> vTrackLength;

   TVector3 VPoint;
   TVector3 VLine1[4];
   TVector3 VLine2[4];
   TGraph2D*  g2D;
   TPolyLine3D *line[4];
   TGraph2D*  g2D_vertex;
   TGraph2D*  g2D_lpoint;
   TH1F * thetaHisto[50];
   TH1F * thetaHisto_coin[50];	
   TH1F * EBeamHisto[50];
   TH1F * ZverHisto[20]; 
   TH1F * thetaCMHisto[50];
   TH1F * ExHisto[50];

};

#endif

#ifdef AnaClass2_cxx
AnaClass2::AnaClass2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("scattering_events_17F_284.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

AnaClass2::~AnaClass2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaClass2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaClass2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaClass2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Event_number = 0;
   Versor_x = 0;
   Versor_y = 0;
   Versor_z = 0;
   Point_x = 0;
   Point_y = 0;
   Point_z = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   Is_primary = 0;
   Total_charge = 0;
   Max_Pos_x = 0;
   Max_Pos_y = 0;
   Max_Pos_z = 0;
   Num_inliers = 0;
   vpx = 0;
   vpy = 0;
   vpz = 0;
   vpt = 0;
   vpe = 0;
   vhits = 0;
   TAC_T = 0;
   TAC_A = 0;
   Trigg_T = 0;
   Trigg_A = 0;
   IC_T = 0;
   IC_A = 0;
   Mesh_T = 0;
   Mesh_A = 0;
   alpha_T = 0;
   alpha_A = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_number", &Event_number, &b_Event_number);
   fChain->SetBranchAddress("Versor_x", &Versor_x, &b_Versor_x);
   fChain->SetBranchAddress("Versor_y", &Versor_y, &b_Versor_y);
   fChain->SetBranchAddress("Versor_z", &Versor_z, &b_Versor_z);
   fChain->SetBranchAddress("Point_x", &Point_x, &b_Point_x);
   fChain->SetBranchAddress("Point_y", &Point_y, &b_Point_y);
   fChain->SetBranchAddress("Point_z", &Point_z, &b_Point_z);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Is_primary", &Is_primary, &b_Is_primary);
   fChain->SetBranchAddress("Total_charge", &Total_charge, &b_Total_charge);
   fChain->SetBranchAddress("Max_Pos_x", &Max_Pos_x, &b_Max_Pos_x);
   fChain->SetBranchAddress("Max_Pos_y", &Max_Pos_y, &b_Max_Pos_y);
   fChain->SetBranchAddress("Max_Pos_z", &Max_Pos_z, &b_Max_Pos_z);
   fChain->SetBranchAddress("Num_inliers", &Num_inliers, &b_Num_inliers);
   fChain->SetBranchAddress("vpx", &vpx, &b_vpx);
   fChain->SetBranchAddress("vpy", &vpy, &b_vpy);
   fChain->SetBranchAddress("vpz", &vpz, &b_vpz);
   fChain->SetBranchAddress("vpt", &vpt, &b_vpt);
   fChain->SetBranchAddress("vpe", &vpe, &b_vpe);
   fChain->SetBranchAddress("vhits", &vhits, &b_vhits);
   fChain->SetBranchAddress("TAC_T", &TAC_T, &b_TAC_T);
   fChain->SetBranchAddress("TAC_A", &TAC_A, &b_TAC_A);
   fChain->SetBranchAddress("Trigg_T", &Trigg_T, &b_Trigg_T);
   fChain->SetBranchAddress("Trigg_A", &Trigg_A, &b_Trigg_A);
   fChain->SetBranchAddress("IC_T", &IC_T, &b_IC_T);
   fChain->SetBranchAddress("IC_A", &IC_A, &b_IC_A);
   fChain->SetBranchAddress("Mesh_T", &Mesh_T, &b_Mesh_T);
   fChain->SetBranchAddress("Mesh_A", &Mesh_A, &b_Mesh_A);
   fChain->SetBranchAddress("alpha_T", &alpha_T, &b_alpha_T);
   fChain->SetBranchAddress("alpha_A", &alpha_A, &b_alpha_A);
   Notify();
}

Bool_t AnaClass2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaClass2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaClass2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

double AnaClass2::sq(double val ){
	
	return val*val;
}




double AnaClass2::omega(double x, double y, double z){
	return sqrt(x*x + y*y + z*z -2*x*y -2*y*z -2*x*z);

}


void AnaClass2::TwobodyKin(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject){

  //in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
  double Et1 = K_proj + m1;
  double Et2 = m2;
  double Et3 = K_eject + m3;
  double Et4  = Et1 + Et2 - Et3;
  double m4_ex;
  double s,t,u; //---Mandelstam variables
  double p1, p3;
  double J_LtoCM; //jacobian Lab to CM

  s = pow(m1,2) + pow(m2,2) +2*m2*Et1;
  u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

  m4_ex = sqrt(  (cos(thetalab) * omega(s,pow(m1,2),pow(m2,2)) * omega(u,pow(m2,2),pow(m3,2)) - (s - pow(m1,2) - pow(m2,2))*(pow(m2,2) + pow(m3,2) - u) )/(2*pow(m2,2)) + s + u - pow(m2,2)  );
  Ex_reco = m4_ex - m4;

  t =   pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;


  //for inverse kinematics Note: this angle corresponds to the recoil
  ThetaCM = 180. - acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) )*180./3.1415 ;

  p1= sqrt(pow(Et1,2)-pow(m1,2));
  p3 = sqrt(pow(Et3,2)-pow(m3,2));

  J_LtoCM = ((omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2)))/(4*s*p1*p3))*(1.+Et1/m2 - cos(thetalab)*(Et3*p1)/(m2*p3)) ;


  //cout<<m1<<"  "<<m2<<"  "<<m3<<endl;
  /*static double output[3];
  output[0]= theta_cm;
  output[1]= Ex;
  output[2]= J_LtoCM;
  return output;
  */
  //return theta_cm;
}



double AnaClass2::getThetaLab(double m1, double m2, double m3, double m4, double K_proj, double thetacm, double Ex){


    double E1 = K_proj + m1;
    double m4_ex = m4 + Ex;
    double s = sq(m1) + sq(m2) +2*m2*E1;
    double a =4.*m2*s;
    double b =(sq(m3)-sq(m4_ex)+s)*(sq(m1)-sq(m2)-s);
    double c = omega(s,sq(m1),sq(m2))*omega(s,sq(m3),sq(m4_ex));

    double E3 = (c*cos(3.141592-thetacm)-b)/a; //estamos viendo el recoil en cm (pi-theta)
    double E4 = (E1 + m2 - E3);
	
    //double K3 = E3 - m3;
    //double K4 = E4 - m4;

    double t =  sq(m2) + sq(m4_ex) - 2*m2*E4;
    double u =  sq(m2) + sq(m3) - 2*m2*E3;



    double thetalab = acos(((s-sq(m1)-sq(m2))*(sq(m2)+sq(m3)-u)+2.*sq(m2)*(sq(m2)+sq(m4_ex)-s-u))/(omega(s,sq(m1),sq(m2))*omega(u,sq(m2),sq(m3))));

	
   return thetalab*180./3.141592;	


}




double AnaClass2::getEkin2(double m1, double m2, double m3, double m4, double K_proj, double thetacm, double Ex){


    double E1 = K_proj + m1;
    double m4_ex = m4 + Ex;
    double s = sq(m1) + sq(m2) +2*m2*E1;
    double a =4.*m2*s;
    double b =(sq(m3)-sq(m4_ex)+s)*(sq(m1)-sq(m2)-s);
    double c = omega(s,sq(m1),sq(m2))*omega(s,sq(m3),sq(m4_ex));

    double E3 = (c*cos(3.141592-thetacm)-b)/a; //estamos viendo el recoil en cm (pi-theta)
    double E4 = (E1 + m2 - E3);
	
    double K3 = E3 - m3;
    //double K4 = E4 - m4;

   
	
   return K3;

}


#endif // #ifdef AnaClass2_cxx
