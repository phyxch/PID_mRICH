#ifndef genMassHypo_h
#define genMassHypo_h

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

using namespace std;

class event;
class hit;
class material;
class Utility;

class genMassHypo
{
 public:
  genMassHypo(string outputfile);
  ~genMassHypo();
  
  int init();
  int process_event(event *aevt, hit *ahit);
  int end();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);
  double probability(TH2D* db, TH2D *hXY);
  
 private:
  material *mat;
  Utility *utility;
  string mOutPutFile;
  TFile *File_OutPut;

  // key: pid | indexSpaceX | indexSpaceY | indexMomentumP | indexMomentumTheta | indexMomentumPhi
  TH1DMap hNEvtvsP; // number of total events
  TH2DMap h_photonDist; // x: photon out_x | y: photon out_y 
};
#endif
