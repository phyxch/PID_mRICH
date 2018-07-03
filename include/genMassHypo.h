#ifndef genMassHypo_h
#define genMassHypo_h

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include "type.h"

using namespace std;

class event;
class hit;
class material;
class Utility;

class genMassHypo
{
 public:
  genMassHypo(string numoflist, string date);
  ~genMassHypo();
  
  int Init();
  int initChain();
  int initHistoMap();
  int initGausSmearing();
  int initHistoQA();

  int Make();

  int Finish();
  int writeHistoMap();
  int writeHistoQA();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);

  double GausSmearing(TF1 *f_gaus);
  bool isInSensorPlane(double out_x, double out_y);
  
 private:
  material *mat;
  Utility *utility;
  string mNumOfList;
  string mDate;
  string mOutPutFile;
  TFile *File_mOutPut;

  TF1 *f_mGaus;
  TH2D *h_mXGausSmearing;
  TH2D *h_mYGausSmearing;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH1DMap h_mNumOfEvents; // number of total events
  TH2DMap h_mPhotonDist; // x: photon out_x | y: photon out_y 
  TH2DMap h_mPhotonGenerated; // x: photon out_x | y: photon out_y 

  TChain *mChainInPut_Events;
  TChain *mChainInPut_Tracks;
};
#endif
