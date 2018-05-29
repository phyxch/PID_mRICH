#ifndef calNSigma_h
#define calNSigma_h

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "type.h"

using namespace std;

class Utility;

class calNSigma
{
 public:
  calNSigma(string date, string outputfile);
  ~calNSigma();
  
  int Init();
  int initHistoMap_Likelihood();
  int initHistoMap_NSigma();

  int Make();

  int Finish();
  int writeHistoMap_NSigma();
  
 private:
  Utility *utility;
  string mDate;
  string mInPutFile;
  string mOutPutFile;
  TFile *mFile_InPut;
  TFile *mFile_OutPut;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH2DMap h_mLikelihoodDiff; // x: p | y: likelihood difference 
  TH1DMap h_mNSigma; // p vs. N_sigma 
};
#endif
