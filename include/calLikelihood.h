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
#include "type.h"

using namespace std;

class event;
class hit;
class material;
class Utility;

class calLikelihood
{
 public:
  calLikelihood(string date, string inputdatabase, string outputfile);
  ~calLikelihood();
  
  int Init();
  int initChain();
  int initTree();
  int initHistoMap();

  int Make();

  int Finish();
  int writeTree();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);
  double probability(TH2D* h_database, TH2D *h_photonDist_PID);
  
 private:
  material *mat;
  Utility *utility;
  string mDate;
  string mInPutDataBase;
  string mOutPutFile;
  TFile *mFile_InPutDataBase;
  TFile *mFile_OutPut;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH1DMap h_mNumOfEvents; // number of total events
  TH2DMap h_mPhotonDist; // x: photon out_x | y: photon out_y 
  
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
  int mNelSbKCs;
  int mNelGaAsP;
  int mNelGaAs;
  double mLpion;
  double mLKaon;
  double mLproton;

  TChain *mChainInPut_Events;
  TChain *mChainInPut_Tracks;
};
#endif
