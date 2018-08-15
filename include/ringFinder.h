#ifndef RingFinder_h
#define RingFinder_h

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "string.h"
#include "type.h"

using namespace std;
class TH1D;
class TH2D;
class TH3D;
class TChain;
class TF1;
class TFile;

class event;
class hit;
class material;
class Utility;

class RingFinder
{
 public:
  RingFinder(string numoflist, string date);
  ~RingFinder();
  
  int Init();
  int initChain();
  int initHistoMap();
  int initGausSmearing();
  int initHistoQA();
  int initHistoCherenkov();

  int Make();
  int clearHistoMap(); 
  int HoughTransform(TH1D *h_NumOfPhotons, TH2D *h_PhotonDist, intVec xPixel, intVec yPixel);
  bool findRing(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit, double &x_Cherenkov, double &y_Cherenkov, double &r_Cherenkov);
  bool isSamePosition(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit);
  bool isCollinear(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit);

  int Finish();
  int writeHistoMap();
  int writeHistoQA();
  int writeHistoCherenkov();
  
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

  // event-by-event analysis
  TH2D *h_mPhotonGenerated; // x: photon out_x | y: photon out_y 
  TH1D *h_mNumOfPhotons; // number of total photons in each event 
  TH2D *h_mPhotonDist; // x: photon out_x | y: photon out_y with detector effect
  intVec mXPixelMap; // corresponding binX number for each photon hit
  intVec mYPixelMap; // corresponding binY number for each photon hit
  TH3D *h_mHoughTransform; // x | y | R
  TH3D *h_mQA_HT; // QA for hough transform

  // all events
  TH3D *h_mCherenkovRing; // x | y | R
  TH1D *h_mNumOfCherenkovPhotons; // 

  TChain *mChainInPut_Events;
  TChain *mChainInPut_Tracks;
};
#endif
