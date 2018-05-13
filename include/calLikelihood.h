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
#include "LLTreeDst.h"
#include "type.h"

const double DEG=180./3.1415926;
const unsigned int nPads=105;   ///// number of photonsenor segmentation pads
const double halfWidth=52.5; // (mm) half width (x) and height (y) of the photon sensor

using namespace std;

class event;
class hit;
class material;
class Utility;

class calLikelihood
{
 public:
  calLikelihood(string inputdatebase, string outputfile);
  ~calLikelihood();
  
  int init();
  int process_event(event *aevt, hit *ahit);
  int end();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);
  double probability(TH2D* h_database, TH2D *h_photonDist_PID);
  
 private:
  material *mat;
  Utility *utility;
  string mInPutDataBase;
  string mOutPutFile;
  TFile *File_InPutDataBase;
  TFile *File_OutPut;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH1DMap hNEvtvsP; // number of total events
  TH2DMap h_photonDist; // x: photon out_x | y: photon out_y 
  
  TTree *mTree;
  LLTreeDst mDst;
};
#endif
