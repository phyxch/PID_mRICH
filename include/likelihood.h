#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "LLTreeDst.h"

const double DEG=180./3.1415926;
//const unsigned int nPads=88;   ///// number of photonsenor segmentation pads
//const double halfWidth=44; // (mm) half width (x) and height (y) of the photon sensor
const unsigned int nPads=8;   ///// number of photonsenor segmentation pads
const double halfWidth=24.25; // (mm) half width (x) and height (y) of the photon sensor

using namespace std;

class event;
class hit;
class material;
class ring;

class likelihood
{
 public:
  likelihood(char *fout);
  ~likelihood();
  
  int init();
  int process_event(event *aevt, hit *ahit);
  int end();
  
  bool isPhoton(hit *ahit, int i);
  bool isReflection(hit *ahit, int i);
  bool isOnAerogel(hit *ahit, int i);
  bool isOnPhotonSensor(hit *ahit, int i);
  void Smearing2D(double inx, double iny, double& outx, double& outy);
  double probability(TH2D* db, TH2D *hXY);
  
 private:
  TRandom *rd;
  material *mat;
  char *outf;
  TFile *outputfile;
  
  TH2D *hNEvtvsP;
  TH3D *hPiXYvsP;
  TH3D *hKaonXYvsP;
  TH3D *hProtonXYvsP;
};
#endif
