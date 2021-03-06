//////////////////////////////////////////////////////////
// RANSAC Class
// J.C. Zamora, jczamorac@gmail.com
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#ifndef Ransac_H
#define Ransac_H

// ROOT Headers
#include <TObject.h>
#include <TH2F.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>

#include <iostream>
#include <chrono>
#include <ctime>
#include <stdio.h>
#include <vector>

using namespace std;



class Ransac
{
  public:
  Ransac();
  ~Ransac();

  void Reset();
	void Init(vector<double> v1, vector<double> v2, vector<double> v3, vector<double> v4);
	void Solve(double dist_thres, double Nminpoints, int Nintera);
  vector<int> RandSam(vector<int> indX);
  void EstimModel(const std::vector<int>  samplesIdx);
  double EstimError(int i);
	vector<double> GetChargeOfTracks();
	vector<double> GetTrackLength();
	double Fit3D(vector<int> inliners, TVector3& V1, TVector3& V2);
  TRandom* Rand;
  TVector3 Vs;
  TVector3 Ps;


  //-----get a cluster (maybe a track class?)
  struct Cluster // return type of structure
    {
      double ClusterStrength;		// strength
      size_t ClusterSize;			// size
      double ClusterChi2;			// size
      std::vector<int> ClusterIndex;			// Indices
      TVector3 ClusterFitP1;			// point 1 from the fitted line
      TVector3 ClusterFitP2;			// point 2 from the fitted line
    };


  typedef std::vector<Cluster> AllClusters;


  void SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2);
  inline AllClusters GetClusters(){return cluster_vector;}
  AllClusters cluster_vector;



  private:
  vector<double> vX, vY, vZ, vQ;
	vector<double> vTrackCharge;
  float fRANSACThreshold;
	float fRANSACPointThreshold;
	float fRANSACChargeThreshold;
	float fRANSACDistance;
	float fRANSACMaxIteration;
	int fNumberOfTracksMax;
	int fOriginalCloudSize;
	double fTotalCharge;
	int fNumberOfPadsX;
	int fNumberOfPadsY;
	int fVerbose;



  public:
  //ClassDef(Ransac,0)
};

#endif
