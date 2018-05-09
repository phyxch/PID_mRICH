#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include "string.h"
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom.h>
#include "TString.h"
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/likelihood.h"

using namespace std;
using namespace TMath;


likelihood::likelihood(char *fout): outf(fout)
{
 cout<<endl;
 cout<<"likelihood::likelihood() ----- Constructor ! ------"<<endl;
 cout<<endl;
 init();
}

likelihood::~likelihood()
{
 cout<<"likelihood::~likelihood() ----- Release memory ! ------"<<endl;
 delete rd;
 delete mat;
 delete outf;
 delete outputfile;
 delete hNEvtvsP;
 delete hPiXYvsP;
 delete hKaonXYvsP;
 delete hProtonXYvsP;
}

int likelihood::init()
{
 cout<<"likelihood::init() ----- Initialization ! ------"<<endl;

 rd = new TRandom();    /// random number generator
 mat = new material(); //// initialize the material

 cout<<"likelihood::init(), create output file: "<< outf <<endl;
 outputfile = new TFile(outf,"RECREATE");

 cout<<"likelihood::init(), initialize database histograms ;"<<endl;
 hNEvtvsP = new TH2D("hNEvtvsP","hNEvtvsP",5,0.,5.,65,2.5,15.5);
 /*
 hPiXYvsP = new TH3D("hPiXYvsP","hPiXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 hKaonXYvsP = new TH3D("hKaonXYvsP","hKaonXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 hProtonXYvsP = new TH3D("hProtonXYvsP","hProtonXYvsP",nPads,-44.,44.,nPads,-44.,44.,65,2.5,15.5);
 */
 hPiXYvsP = new TH3D("hPiXYvsP","hPiXYvsP",105,-52.5,52.5,105,-52.5,52.5,65,2.5,15.5);
 hKaonXYvsP = new TH3D("hKaonXYvsP","hKaonXYvsP",105,-52.5,52.5,105,-52.5,52.5,65,2.5,15.5);
 hProtonXYvsP = new TH3D("hProtonXYvsP","hProtonXYvsP",105,-52.5,52.5,105,-52.5,52.5,65,2.5,15.5);
 
 return 0;
}

int likelihood::process_event(event *aevt, hit *ahit)
{
  int pid_gen=0;
  float px_gen=0;
  float py_gen=0;
  float pz_gen=0;
  float vx_gen=0;
  float vy_gen=0;
  float vz_gen=0;
  float p_gen=0;
  float theta_gen=0;
  float phi_gen=0;
  // cout<< "tree generated size: "<< aevt->get_pid()->size() <<";    tree flux size:  "<< ahit->get_hitn()->size() <<endl;
  for (unsigned int i=0;i<aevt->get_pid()->size();i++) {
    // cout << aevt->get_pid()->at(i) << " " << aevt->get_px()->at(i) << " " << aevt->get_py()->at(i) << " " << aevt->get_pz()->at(i) << " " << aevt->get_vx()->at(i) << " " << aevt->get_vy()->at(i) << " " << aevt->get_vz()->at(i) << endl; 
    pid_gen=aevt->get_pid()->at(i);
    px_gen=aevt->get_px()->at(i)/1e3;    //in MeV, convert to GeV
    py_gen=aevt->get_py()->at(i)/1e3;    //in MeV, convert to GeV
    pz_gen=aevt->get_pz()->at(i)/1e3;    //in MeV, convert to GeV
    vx_gen=aevt->get_vx()->at(i)/1e1;    //in mm, convert to cm
    vy_gen=aevt->get_vy()->at(i)/1e1;    //in mm, convert to cm
    vz_gen=aevt->get_vz()->at(i)/1e1;    //in mm, convert to cm
    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    theta_gen=acos(pz_gen/p_gen)*DEG;    //in deg
    phi_gen=atan2(py_gen,px_gen)*DEG;    //in deg            
  }  
  
  if(pid_gen==-211) hNEvtvsP->Fill(0.,p_gen);
  else if(pid_gen==-321) hNEvtvsP->Fill(1.,p_gen);
  else if(pid_gen==2212) hNEvtvsP->Fill(2.,p_gen);
  
  int nhits = ahit->get_hitn()->size();
  for (int i=0;i<nhits;i++) {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i)){
      double out_x(0.), out_y(0.);
      Smearing2D(ahit->get_out_x()->at(i), ahit->get_out_y()->at(i), out_x, out_y);
      if(fabs(out_x)>halfWidth||fabs(out_y)>halfWidth) continue;
      double photonE = ahit->get_trackE()->at(i);   /// in MeV (GEANT4 default)
      double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
      double QE_GaAsP=mat->extrapQE_GaAsP(wavelength);
      if(QE_GaAsP>rd->Uniform(0.,1.)){
	if(pid_gen==-211) hPiXYvsP->Fill(out_x,out_y,p_gen);
	else if(pid_gen==-321) hKaonXYvsP->Fill(out_x,out_y,p_gen);
	else if(pid_gen==2212) hProtonXYvsP->Fill(out_x,out_y,p_gen);
      }
    }
  }
  
  return 0;
}

int likelihood::end()
{
  cout<<endl;
  cout<<"likelihood::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(outputfile != NULL){
    outputfile->cd();
    hNEvtvsP->Write();
    hPiXYvsP->Write();
    hKaonXYvsP->Write();
    hProtonXYvsP->Write();
    outputfile->Close();
  }
  return 0;
}

bool likelihood::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool likelihood::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool likelihood::isOnAerogel(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=55.5 && ahit->get_out_z()->at(i)<=85.5) return true;
  //if(ahit->get_out_z()->at(i)>=50.5 && ahit->get_out_z()->at(i)<=70.5) return true;
  else return false;
}

bool likelihood::isOnPhotonSensor(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=159.5) return true;
  //if(ahit->get_out_z()->at(i)>143 && ahit->get_out_z()->at(i)<147) return true;
  else return false;
}

double likelihood::probability(TH2D *db, TH2D *hXY)
{
  double prob=1.0;
  const double noise = 2./nPads/nPads; ///// noise level is set to be 2 electrons on photonsensor
  
  for(unsigned int i=1; i<=nPads; i++){   ///// photonsensor pad loop in X
    for(unsigned int j=1; j<=nPads; j++)   ///// photonsensot pad loop in Y
      {
	double k = hXY->GetBinContent(i,j); ///// should be integer (mostly 0 or 1)
	double lambda = db->GetBinContent(i,j); ///// get the expected value of Poisson distribution from DB
	if(lambda==0.) prob *= PoissonI(k,noise);
	else if(lambda>0.) prob *= PoissonI(k,lambda);
	else cout<<"DB histogram bin content should not be < 0. !"<<endl;
      }
  }
  return log(prob);
}

void likelihood::Smearing2D(double inx, double iny, double& outx, double& outy)
{
  TF1 *fx = new TF1("fx","1/([1]*sqrt(2*3.1415926))*exp(-1*pow(x-[0],2)/(2.*[1]*[1]))",inx-20.,inx+20.);
  fx->SetParameter(0,inx);
  fx->SetParameter(1,1.); //// 1 mm smearing in photon position x
  outx = fx->GetRandom();
  
  TF1 *fy = new TF1("fy","1/([1]*sqrt(2*3.1415926))*exp(-1*pow(x-[0],2)/(2.*[1]*[1]))",iny-20.,iny+20.);
  fy->SetParameter(0,iny);
  fy->SetParameter(1,1.); //// 1 mm smearing in photon position y
  outy = fy->GetRandom();
  
  delete fx;
  delete fy;
  return ;
}


////// This is the main function 
int main(int argc, char **argv)
{
  TString inputdir = "/work/eic/xusun/output/modular_rich/May02_2018/";

  TString InPutList = "/work/eic/xusun/list/modular_rich/mRICH_test.list";

  char outputf[256];
  strcpy(outputf,"/work/eic/xusun/output/modular_rich/test/test.root");
  
  TChain *fevt  = new TChain("generated");
  TChain *fhit  = new TChain("eic_rich");

  if (!InPutList.IsNull())   // if input file is ok
  {
    TString InFo_List ="Open test file list ";
    cout << InFo_List.Data() << endl;
    ifstream in(InPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = inputdir+addfile;
	  fevt->AddFile(addfile.Data(),-1,"generated");
	  fhit->AddFile(addfile.Data(),-1,"eic_rich");
	  Long64_t file_entries = fevt->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      TString InFo_Warning = "WARNING: test file input is problemtic";
      cout << InFo_Warning.Data() << endl;
    }
  }
  
  event *aevt = new event(fevt);  /// declear and save info to branchs for event
  hit *ahit = new hit(fhit);	/// declear and save info to branchs for track
  
  likelihood *lk = new likelihood(outputf);
  
  int nevent = (int)fevt->GetEntries();
  cout << "total number of events:  " << nevent << endl;
  for (Int_t i=0;i<nevent;i++) { 
    if(i%100==0) cout << "processing event:  " << i << " ;"<<endl;
    fevt->GetEntry(i);  
    fhit->GetEntry(i);
    lk->process_event(aevt, ahit);
  }
  
  lk->end();
  delete ahit;
  delete aevt;
  delete fhit;
  delete fevt;
  
  return 0;
}
