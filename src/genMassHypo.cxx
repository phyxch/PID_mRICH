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
#include "../include/genMassHypo.h"

using namespace std;
using namespace TMath;


genMassHypo::genMassHypo(char *fout): outf(fout)
{
 cout<<endl;
 cout<<"genMassHypo::genMassHypo() ----- Constructor ! ------"<<endl;
 cout<<endl;
 init();
}

genMassHypo::~genMassHypo()
{
 cout<<"genMassHypo::~genMassHypo() ----- Release memory ! ------"<<endl;
 delete mat;
 delete outf;
 delete outputfile;
 delete hNEvtvsP;
 delete h_photonDist_piplus;
 delete h_photonDist_piminus;
 delete h_photonDist_Kplus;
 delete h_photonDist_Kminus;
 delete h_photonDist_proton;
 delete h_photonDist_antiproton;
}

int genMassHypo::init()
{
 cout<<"genMassHypo::init() ----- Initialization ! ------"<<endl;

 mat = new material(); //// initialize the material

 cout<<"genMassHypo::init(), create output file: "<< outf <<endl;
 outputfile = new TFile(outf,"RECREATE");

 cout<<"genMassHypo::init(), initialize database histograms ;"<<endl;
 hNEvtvsP = new TH2D("hNEvtvsP","hNEvtvsP",6,-3.0,3.0,65,2.5,15.5);

 h_photonDist_piplus = new TH3D("h_photonDist_piplus","h_photonDist_piplus",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 h_photonDist_piminus = new TH3D("h_photonDist_piminus","h_photonDist_piminus",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);

 h_photonDist_Kplus = new TH3D("h_photonDist_Kplus","h_photonDist_Kplus",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 h_photonDist_Kminus = new TH3D("h_photonDist_Kminus","h_photonDist_Kminus",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 
 h_photonDist_proton = new TH3D("h_photonDist_proton","h_photonDist_proton",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);
 h_photonDist_antiproton = new TH3D("h_photonDist_antiproton","h_photonDist_antiproton",nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth,65,2.5,15.5);

 return 0;
}

int genMassHypo::process_event(event *aevt, hit *ahit)
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
  
  if(pid_gen==-211) hNEvtvsP->Fill(-2.5,p_gen);
  else if(pid_gen==-321) hNEvtvsP->Fill(-1.5,p_gen);
  else if(pid_gen==-2212) hNEvtvsP->Fill(-0.5,p_gen);
  if(pid_gen==211) hNEvtvsP->Fill(2.5,p_gen);
  else if(pid_gen==321) hNEvtvsP->Fill(1.5,p_gen);
  else if(pid_gen==2212) hNEvtvsP->Fill(0.5,p_gen);
  
  int nhits = ahit->get_hitn()->size();
  for (int i=0;i<nhits;i++) 
  {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i))
    {
      double out_x = ahit->get_out_x()->at(i);
      double out_y = ahit->get_out_y()->at(i);
      if(pid_gen==211)        h_photonDist_piplus->Fill(out_x,out_y,p_gen);
      else if(pid_gen==-211)  h_photonDist_piminus->Fill(out_x,out_y,p_gen);
      else if(pid_gen==321)   h_photonDist_Kplus->Fill(out_x,out_y,p_gen);
      else if(pid_gen==-321)  h_photonDist_Kminus->Fill(out_x,out_y,p_gen);
      else if(pid_gen==2212)  h_photonDist_proton->Fill(out_x,out_y,p_gen);
      else if(pid_gen==-2212) h_photonDist_antiproton->Fill(out_x,out_y,p_gen);
    }
  }
  
  return 0;
}

int genMassHypo::end()
{
  cout<<endl;
  cout<<"genMassHypo::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(outputfile != NULL){
    outputfile->cd();
    hNEvtvsP->Write();
    h_photonDist_piplus->Write();
    h_photonDist_piminus->Write();
    h_photonDist_Kplus->Write();
    h_photonDist_Kminus->Write();
    h_photonDist_proton->Write();
    h_photonDist_antiproton->Write();
    outputfile->Close();
  }
  return 0;
}

bool genMassHypo::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool genMassHypo::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool genMassHypo::isOnAerogel(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=55.5 && ahit->get_out_z()->at(i)<=85.5) return true;
  //if(ahit->get_out_z()->at(i)>=50.5 && ahit->get_out_z()->at(i)<=70.5) return true;
  else return false;
}

bool genMassHypo::isOnPhotonSensor(hit *ahit, int i)
{
  double out_z = ahit->get_out_z()->at(i);
  if(out_z > 253.6624 && out_z < 255.1626) return true;
  else return false;
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
  
  genMassHypo *genMassHypotheses = new genMassHypo(outputf);
  
  int nevent = (int)fevt->GetEntries();
  cout << "total number of events:  " << nevent << endl;
  for (Int_t i=0;i<nevent;i++) { 
    if(i%100==0) cout << "processing event:  " << i << " ;"<<endl;
    fevt->GetEntry(i);  
    fhit->GetEntry(i);
    genMassHypotheses->process_event(aevt, ahit);
  }
  
  genMassHypotheses->end();
  delete ahit;
  delete aevt;
  delete fhit;
  delete fevt;
  
  return 0;
}
