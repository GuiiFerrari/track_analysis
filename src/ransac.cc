//////////////////////////////////////////////////////////
// RANSAC Class
// J.C. Zamora, jczamorac@gmail.com
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#define RANSAC_cc

#include "ransac.h"
#include<iostream>

using namespace std;

//ClassImp(Ransac)

Ransac::Ransac() {

  fRANSACMaxIteration = 500;
	fRANSACThreshold = 30;//30;  //min number of points
	fRANSACPointThreshold = 0.10;//0.07;
	fRANSACChargeThreshold = 0.10;//0.07;
	fRANSACDistance = 15;//7;
}

Ransac::~Ransac() {

  delete Rand;
}

void Ransac::Init(vector<double> v1, vector<double> v2, vector<double> v3, vector<double> v4)
{

  Reset();
  vX = v1;
	vY = v2;
	vZ = v3;
	vQ = v4;

	fOriginalCloudSize = vX.size();
	double TotalCharge=0;
	for(unsigned int i=0; i< vQ.size(); i++){
		TotalCharge += vQ[i];
	}
	fTotalCharge = TotalCharge;

  Rand =new TRandom();

  //std::cout << "/* Ransac inicializa  */" << '\n';

}

void Ransac::Reset()
{

	vX.clear();
	vY.clear();
	vZ.clear();
	vQ.clear();
  Vs.Clear();
  Ps.Clear();
  //delete Rand;

}


void Ransac::Solve(double dist_thres, double Nminpoints, int Nintera)
{

    fRANSACDistance = dist_thres;
    fRANSACThreshold = Nminpoints;
    fRANSACMaxIteration = Nintera;

  	double RemainingCharge = fTotalCharge;
  	//vector<double> trackX, trackY, trackZ;
  	//vector<double> trackY;

    std::vector<int> remainIndex;
    for (size_t i = 0; i < vX.size(); i++)
    remainIndex.push_back(i);



  	  double minErr=1E20;
  		//double track_charge=0;
  		//trackX.clear();
  		//trackY.clear();
  		//trackZ.clear();

  	  TVector3 V1, V2;
  	  std::vector< int> inliners;
  		inliners.clear();



  	  for(int i=0;i<fRANSACMaxIteration;i++){

        if(remainIndex.size()<fRANSACThreshold) break;

        std::vector< int> Rsamples = RandSam(remainIndex);
        EstimModel(Rsamples);

  	    //std::vector<int> alsoin;
  	   	//alsoin.clear();
        std::vector<int> inlIdxR;
        int nbInliers = 0;
        double weight = 0;

        for (auto j = remainIndex.begin(); j != remainIndex.end(); ++j){

  	      double error = EstimError(*j);
          error = error*error;
  	      if(error<(fRANSACDistance*fRANSACDistance)){
  	        	inlIdxR.push_back(*j);
              nbInliers++;
              weight +=error;
  	        	}

  		}


  	    if(nbInliers>fRANSACThreshold){
  	    	TVector3 v1, v2;
  	      double chi2=Fit3D(inlIdxR,v1,v2);
          SetCluster(inlIdxR, weight/nbInliers, chi2,v1,v2);
  	      //if(chi2<minErr){
  	      	//minErr=chi2;
  	        //V1=v1;
  	        //V2=v2;
  	        //inliners=alsoin;

            std::vector<int> tempRemain;
            std::set_difference(remainIndex.begin(), remainIndex.end(), inlIdxR.begin(), inlIdxR.end(),
            std::inserter(tempRemain, tempRemain.begin()));
            remainIndex = tempRemain;
            tempRemain.clear();
            inlIdxR.clear();
  	     	//}
  			}

  		}//for RANSAC interactions
}


vector<int> Ransac::RandSam(vector<int> indX)
{
  size_t pclouds = indX.size();
  int p1,p2;
  /*
  //-------Uniform sampling
  p1=(int)(Rand->Uniform(0,pclouds));

	 do{
      p2=(int)(Rand->Uniform(0,pclouds));
    } while(p2==p1);

  vector<int> ranpair{indX[p1],indX[p2]};
  */

  //--------Gaussian sampling
  double dist = 0;
  double sigma = 10.0;
  double y = 0;
  double gauss = 0;
	 do{
      p1=(int)(Rand->Uniform(0,pclouds));
      p2=(int)(Rand->Uniform(0,pclouds));
      TVector3 P1 ={vX[indX[p1]],vY[indX[p1]],vZ[indX[p1]]};
      TVector3 P2 ={vX[indX[p2]],vY[indX[p2]],vZ[indX[p2]]};
      TVector3 dif = P2-P1;
      dist = dif.Mag();
      gauss = 1.0/(sigma*sqrt(2.0*TMath::Pi()))*exp(-0.5*pow(dist/sigma,2));
      y = (Rand->Uniform(0,1));
    } while(p2==p1 || y>gauss);

  vector<int> ranpair{indX[p1],indX[p2]};



  return ranpair;

}


void Ransac::EstimModel(const std::vector<int>  samplesIdx)
{

  //line from two points
  TVector3 Po1 = {vX[samplesIdx[0]], vY[samplesIdx[0]], vZ[samplesIdx[0]]};
  TVector3 Po2 = {vX[samplesIdx[1]], vY[samplesIdx[1]], vZ[samplesIdx[1]]};

  Vs = Po2 - Po1;
  Ps = Po1;

}

double Ransac::EstimError(int i)
{
    //distance point to line
    TVector3 newPoint = {vX[i], vY[i], vZ[i]};
    TVector3 vec = Ps - newPoint;
    TVector3 nD = Vs.Cross(vec);

	  double dist = nD.Mag()/Vs.Mag();


    return  dist;
}


void Ransac::SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2)
{

    Cluster cstr;
    cstr.ClusterIndex = samplesIdx;
    cstr.ClusterSize = samplesIdx.size();
    cstr.ClusterStrength = cost;
    cstr.ClusterChi2 = Chi2;
    cstr.ClusterFitP1 = CP1;
    cstr.ClusterFitP2 = CP2;
    cluster_vector.push_back(cstr);
}


double Ransac::Fit3D( vector<int> inliners, TVector3& V1, TVector3& V2)
{
    //------3D Line Regression
    //----- adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D
    int R, C;
    double Q;
    double Xm,Ym,Zm;
    double Xh,Yh,Zh;
    double a,b;
    double Sxx,Sxy,Syy,Sxz,Szz,Syz;
    double theta;
    double K11,K22,K12,K10,K01,K00;
    double c0,c1,c2;
    double p,q,r,dm2;
    double rho,phi;

    Q=Xm=Ym=Zm=0.;
		double total_charge=0;
    Sxx=Syy=Szz=Sxy=Sxz=Syz=0.;

    for (auto i : inliners)
    {
        Q+=vQ[i]/10.;
        Xm+=vX[i]*vQ[i]/10.;
        Ym+=vY[i]*vQ[i]/10.;
        Zm+=vZ[i]*vQ[i]/10.;
        Sxx+=vX[i]*vX[i]*vQ[i]/10.;
        Syy+=vY[i]*vY[i]*vQ[i]/10.;
        Szz+=vZ[i]*vZ[i]*vQ[i]/10.;
        Sxy+=vX[i]*vY[i]*vQ[i]/10.;
        Sxz+=vX[i]*vZ[i]*vQ[i]/10.;
        Syz+=vY[i]*vZ[i]*vQ[i]/10.;
    }
    //vTrackCharge.push_back(total_charge);

    Xm/=Q;
    Ym/=Q;
    Zm/=Q;
    Sxx/=Q;
    Syy/=Q;
    Szz/=Q;
    Sxy/=Q;
    Sxz/=Q;
    Syz/=Q;
    Sxx-=(Xm*Xm);
    Syy-=(Ym*Ym);
    Szz-=(Zm*Zm);
    Sxy-=(Xm*Ym);
    Sxz-=(Xm*Zm);
    Syz-=(Ym*Zm);

    theta=0.5*atan((2.*Sxy)/(Sxx-Syy));

    K11=(Syy+Szz)*pow(cos(theta),2)+(Sxx+Szz)*pow(sin(theta),2)-2.*Sxy*cos(theta)*sin(theta);
    K22=(Syy+Szz)*pow(sin(theta),2)+(Sxx+Szz)*pow(cos(theta),2)+2.*Sxy*cos(theta)*sin(theta);
    K12=-Sxy*(pow(cos(theta),2)-pow(sin(theta),2))+(Sxx-Syy)*cos(theta)*sin(theta);
    K10=Sxz*cos(theta)+Syz*sin(theta);
    K01=-Sxz*sin(theta)+Syz*cos(theta);
    K00=Sxx+Syy;

    c2=-K00-K11-K22;
    c1=K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
    c0=K01*K01*K11+K10*K10*K22-K00*K11*K22;


    p=c1-pow(c2,2)/3.;
    q=2.*pow(c2,3)/27.-c1*c2/3.+c0;
    r=pow(q/2.,2)+pow(p,3)/27.;


    if(r>0) dm2=-c2/3.+pow(-q/2.+sqrt(r),1./3.)+pow(-q/2.-sqrt(r),1./3.);
    if(r<0)
    {
        rho=sqrt(-pow(p,3)/27.);
        phi=acos(-q/(2.*rho));
        dm2=min(-c2/3.+2.*pow(rho,1./3.)*cos(phi/3.),min(-c2/3.+2.*pow(rho,1./3.)*cos((phi+2.*TMath::Pi())/3.),-c2/3.+2.*pow(rho,1./3.)*cos((phi+4.*TMath::Pi())/3.)));
    }

    a=-K10*cos(theta)/(K11-dm2)+K01*sin(theta)/(K22-dm2);
    b=-K10*sin(theta)/(K11-dm2)-K01*cos(theta)/(K22-dm2);

    Xh=((1.+b*b)*Xm-a*b*Ym+a*Zm)/(1.+a*a+b*b);
    Yh=((1.+a*a)*Ym-a*b*Xm+b*Zm)/(1.+a*a+b*b);
    Zh=((a*a+b*b)*Zm+a*Xm+b*Ym)/(1.+a*a+b*b);

    V1.SetXYZ(Xm,Ym,Zm);
    V2.SetXYZ(Xh,Yh,Zh);

    return(fabs(dm2/Q));
}
