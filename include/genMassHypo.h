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

const double DEG=180./3.1415926;
//const unsigned int nPads=88;   ///// number of photonsenor segmentation pads
//const double halfWidth=44; // (mm) half width (x) and height (y) of the photon sensor
const unsigned int nPads=105;   ///// number of photonsenor segmentation pads
const double halfWidth=52.5; // (mm) half width (x) and height (y) of the photon sensor

using namespace std;

class event;
class hit;
class material;

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
  void Smearing2D(double inx, double iny, double& outx, double& outy);
  double probability(TH2D* db, TH2D *hXY);
  
 private:
  material *mat;
  string mOutPutFile;
  TFile *File_OutPut;
  
  TH2D *hNEvtvsP;
  TH3D *h_photonDist_piplus; // x: photon out_x | y: photon out_y | z: total momentum of generated particle
  TH3D *h_photonDist_piminus;
  TH3D *h_photonDist_Kplus;
  TH3D *h_photonDist_Kminus;
  TH3D *h_photonDist_proton;
  TH3D *h_photonDist_antiproton;
};
#endif