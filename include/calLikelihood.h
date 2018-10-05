#ifndef calLikelihood_h
#define calLikelihood_h

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include "type.h"

using namespace std;

class event;
class hit;
class material;
class Utility;

class calLikelihood
{
 public:
  calLikelihood(string numoflist, string date, string inputdatabase);
  ~calLikelihood();
  
  int Init();
  int initChain();
  int initTree();
  int initHistoMap();
  int initHistoQA();
  int initGausSmearing();

  int Make();

  int Finish();
  int writeTree();
  int writeQA();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);

  double probability(TH2D* h_database, TH2D *h_photonDist_PID);
  double GausSmearing(TF1 *f_gaus);
  bool isInSensorPlane(double out_x, double out_y);
  
 private:
  material *mat;
  Utility *utility;
  string mNumOfList;
  string mDate;
  string mInPutDataBase;
  string mOutPutFile;
  TFile *mFile_InPutDataBase;
  TFile *mFile_OutPut;

  TF1 *f_mGaus;
  TH2D *h_mXGausSmearing;
  TH2D *h_mYGausSmearing;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH1DMap h_mNumOfEvents; // number of total events
  TH2DMap h_mPhotonDist; // x: photon out_x | y: photon out_y 

  TH2D *h_mPhotonEnergy_Wavelength;
  TH2D *hnHitAeglPerEvtvsMom;
  TH2D *hnPhotonElvsnHits_SbKCs;
  TH2D *hnPhotonElvsnHits_GaAsP;
  TH2D *hnPhotonElvsnHits_GaAs;
  
  TTree *mTree;
  int mPid;
  double mPx;
  double mPy;
  double mPz;
  double mVx;
  double mVy;
  double mVz;
  double mTheta;
  double mPhi;
  int mNHit;
  int mNHitAegl;
  int mNHitPhoDet;
  int mNElSbKCs;
  int mNElGaAsP;
  int mNElGaAs;
  double mLpion;
  double mLKaon;
  double mLproton;
  double mLelectron;

  TChain *mChainInPut_Events;
  TChain *mChainInPut_Tracks;
};
#endif
