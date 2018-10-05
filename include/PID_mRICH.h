#ifndef PID_mRICH_h
#define PID_mRICH_h

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include <utility>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "type.h"

using namespace std;

class Utility;

class PID_mRICH
{
 public:
  PID_mRICH(string date, string outputfile);
  ~PID_mRICH();
  
  int Init();
  int initChain();
  int initHistoMap_Likelihood();
  int initHistoMap_Probability();

  int Make();
  std::pair<double,double> get_LikelihoodDiff(int pid, double Lpion, double LKaon, double Lproton, double Lelectron);
  std::pair<int,int> get_misPID(int pid);
  int get_rank(int pid, double Lpion, double LKaon, double Lproton, double Lelectron);

  int Finish();
  int writeHistoMap_Likelihood();
  int writeHistoMap_Probability();
  
 private:
  string mDate;
  string mOutPutFile;
  TFile *mFile_OutPut;
  
  Utility *utility;

  TChain *mChainInPut; // branches list below
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
  double mLelectron;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumTheta | indexMomentumPhi
  TH1DMap h_mProbability; // p vs. probability
  TH1DMap h_mSumOfPID; // p vs. total number
  TH2DMap h_mLikelihoodDiff; // x: p | y: likelihood difference 
  // doubleMap SumOfPID;
  
};
#endif
