#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include <utility>
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
#include "../include/Utility.h"

using namespace std;
using namespace TMath;


genMassHypo::genMassHypo(string outputfile)
{
 cout<<endl;
 cout<<"genMassHypo::genMassHypo() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mOutPutFile = outputfile;
 init();
}

genMassHypo::~genMassHypo()
{
 cout<<"genMassHypo::~genMassHypo() ----- Release memory ! ------"<<endl;
 delete mat;
 delete utility;
 delete File_OutPut;
}

int genMassHypo::init()
{
 cout<<"genMassHypo::init() ----- Initialization ! ------"<<endl;

 mat = new material(); //// initialize the material
 utility = new Utility(); // initialize utility class

 cout<<"genMassHypo::init(), create output file: "<< mOutPutFile.c_str() <<endl;
 File_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

 cout<<"genMassHypo::init(), initialize database histograms ;"<<endl;

 std::string PID[6] = {"piplus","Kplus","proton","piminus","Kminus","antiproton"};
 for(int i_pid = 0; i_pid < 6; ++i_pid)
 {
   for(int i_vx = 0; i_vx < 5; ++i_vx)
   {
     for(int i_vy = 0; i_vy < 5; ++i_vy)
     {
       for(int i_mom = 0; i_mom < 10; ++i_mom)
       {
	 string key_events = Form("h_NumofEvents_%s_vx_%d_vy_%d_mom_%d",PID[i_pid].c_str(),i_vx,i_vy,i_mom);
	 // cout << "genMassHypo::init(), initialize histogram: " << key_events.c_str() << endl;
	 hNEvtvsP[key_events] = new TH1D(key_events.c_str(),key_events.c_str(),1,-0.5,0.5);

	 string key_photon = Form("h_photonDist_%s_vx_%d_vy_%d_mom_%d",PID[i_pid].c_str(),i_vx,i_vy,i_mom);
	 // cout << "genMassHypo::init(), initialize histogram: " << key_photon.c_str() << endl;
	 h_photonDist[key_photon] = new TH2D(key_photon.c_str(),key_photon.c_str(),nPads,-halfWidth,halfWidth,nPads,-halfWidth,halfWidth);
       }
     }
   }
 }

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
    vx_gen=aevt->get_vx()->at(i);        //in mm
    vy_gen=aevt->get_vy()->at(i);        //in mm
    vz_gen=aevt->get_vz()->at(i);        //in mm
    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    theta_gen=acos(pz_gen/p_gen)*DEG;    //in deg
    phi_gen=atan2(py_gen,px_gen)*DEG;    //in deg            
  }  
  
  string identifiedParticle = utility->get_IdentifiedParticle(pid_gen);
  int indexSpaceX = utility->get_indexSpaceX(vx_gen);
  int indexSpaceY = utility->get_indexSpaceX(vy_gen);
  int indexMomentumP = utility->get_indexMomentumP(px_gen,py_gen,pz_gen);

  string key_events = Form("h_NumofEvents_%s_vx_%d_vy_%d_mom_%d",identifiedParticle.c_str(),indexSpaceX,indexSpaceY,indexMomentumP);
  // cout << "fill histogram: " << key_events.c_str() << endl;
  hNEvtvsP[key_events]->Fill(0);

  string key_photon = Form("h_photonDist_%s_vx_%d_vy_%d_mom_%d",identifiedParticle.c_str(),indexSpaceX,indexSpaceY,indexMomentumP);
  // cout << "fill histogram: " << key_photon.c_str() << endl;
  int nhits = ahit->get_hitn()->size();
  for (int i=0;i<nhits;i++) 
  {
    if(isPhoton(ahit,i) && !isReflection(ahit,i) && isOnPhotonSensor(ahit,i))
    {
      double out_x = ahit->get_out_x()->at(i);
      double out_y = ahit->get_out_y()->at(i);
      h_photonDist[key_photon]->Fill(out_x,out_y);
    }
  }

  return 0;
}

int genMassHypo::end()
{
  cout<<endl;
  cout<<"genMassHypo::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(File_OutPut != NULL){
    File_OutPut->cd();
    std::string PID[6] = {"piplus","Kplus","proton","piminus","Kminus","antiproton"};
    for(int i_pid = 0; i_pid < 6; ++i_pid)
    {
      for(int i_vx = 0; i_vx < 5; ++i_vx)
      {
	for(int i_vy = 0; i_vy < 5; ++i_vy)
	{
	  for(int i_mom = 0; i_mom < 10; ++i_mom)
	  {
	    string key_events = Form("h_NumofEvents_%s_vx_%d_vy_%d_mom_%d",PID[i_pid].c_str(),i_vx,i_vy,i_mom);
	    // cout << "genMassHypo::end(), write histogram: " << key_events.c_str() << endl;
	    hNEvtvsP[key_events]->Write();

	    string key_photon = Form("h_photonDist_%s_vx_%d_vy_%d_mom_%d",PID[i_pid].c_str(),i_vx,i_vy,i_mom);
	    // cout << "genMassHypo::end(), write histogram: " << key_photon.c_str() << endl;
	    h_photonDist[key_photon]->Write();
	  }
	}
      }
    }
    File_OutPut->Close();
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
  string date = "May10_2018";

  string inputdir = Form("/work/eic/xusun/output/modular_rich/%s/",date.c_str());
  string InPutList = Form("/work/eic/xusun/list/modular_rich/mRICH_PDF_%s.list",date.c_str());

  string outputfile = Form("/work/eic/xusun/output/database/PDF_database_%s.root",date.c_str());
  
  TChain *fevt  = new TChain("generated");
  TChain *fhit  = new TChain("eic_rich");

  if (!InPutList.empty())   // if input file is ok
  {
    TString InFo_List ="Open test file list ";
    cout << InFo_List.Data() << endl;
    ifstream in(InPutList.c_str());  // input stream
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
  
  genMassHypo *genMassHypotheses = new genMassHypo(outputfile);
  
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
