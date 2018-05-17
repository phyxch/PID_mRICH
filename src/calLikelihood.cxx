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
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/calLikelihood.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;

calLikelihood::calLikelihood(string date, string inputdatabase, string outputfile)
{
 cout<<endl;
 cout<<"calLikelihood::calLikelihood() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mInPutDataBase = inputdatabase;
 mOutPutFile = outputfile;
 mDate = date;
 utility = new Utility(); // initialize utility class
}

calLikelihood::~calLikelihood()
{
 cout<<"calLikelihood::~calLikelihood() ----- Release memory ! ------"<<endl;
 delete mat;
 delete utility;
 delete mFile_InPutDataBase;
 delete mFile_OutPut;
 delete mTree;
}

int calLikelihood::Init()
{
  cout<<"calLikelihood::Init() ----- Initialization ! ------"<<endl;

  cout<<"calLikelihood::Init(), create output file: "<< mOutPutFile.c_str() <<endl;
  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
  mat = new material(); //// initialize the material

  initChain();
  initHistoMap();
  initTree();
  return 0;
}

int calLikelihood::initChain()
{
  string inputdir = Form("/work/eic/xusun/output/modular_rich/%s/",mDate.c_str());
  string InPutList = Form("/work/eic/xusun/list/modular_rich/mRICH_PDF_%s.list",mDate.c_str());
  
  mChainInPut_Events = new TChain("generated");
  mChainInPut_Tracks = new TChain("eic_rich");

  if (!InPutList.empty())   // if input file is ok
  {
    cout << "Open test file list " << endl;
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
	  string addfile;
	  addfile = str;
	  addfile = inputdir+addfile;
	  mChainInPut_Events->AddFile(addfile.c_str(),-1,"generated");
	  mChainInPut_Tracks->AddFile(addfile.c_str(),-1,"eic_rich");
	  long file_entries = mChainInPut_Events->GetEntries();
	  cout << "File added to data chain: " << addfile.c_str() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: test file input is problemtic" << endl;
    }
  }
  return 0;
}

int calLikelihood::initHistoMap()
{
  cout<<"calLikelihood::init() ----- Initialization ! ------"<<endl;

  cout<<"calLikelihood::init(), read database file: "<< mInPutDataBase.c_str() <<endl;
  mFile_InPutDataBase = TFile::Open(mInPutDataBase.c_str());
  for(int i_pid = 0; i_pid < 6; ++i_pid)
  {
    for(int i_vx = 0; i_vx < mRICH::mNumOfIndexSpaceX; ++i_vx)
    {
      for(int i_vy = 0; i_vy < mRICH::mNumOfIndexSpaceY; ++i_vy)
      {
	for(int i_mom = 0; i_mom < mRICH::mNumOfIndexMomentumP; ++i_mom)
	{
	  for(int i_theta = 0; i_theta < mRICH::mNumOfIndexMomentumTheta; ++i_theta)
	  {
	    for(int i_phi = 0; i_phi < mRICH::mNumOfIndexMomentumPhi; ++i_phi)
	    {
	      string key_events = utility->gen_KeyNumOfEvents(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "calLikelihood::init(), initialize histogram: " << key_events.c_str() << endl;
	      h_mNumOfEvents[key_events] = (TH1D*)mFile_InPutDataBase->Get(key_events.c_str())->Clone();

	      string key_photon = utility->gen_KeyMassHypo(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "calLikelihood::init(), initialize histogram: " << key_photon.c_str() << endl;
	      h_mPhotonDist[key_photon] = (TH2D*)mFile_InPutDataBase->Get(key_photon.c_str())->Clone();
	    }
	  }
	}
      }
    }
  }
}

int calLikelihood::initTree()
{
  cout<<"calLikelihood::init(), initialize tree  ; "<<endl;
  mTree = new TTree("LLTreeDst","Tree for Likelihood Analysis");
  mTree->SetDirectory(mFile_OutPut);
  mTree->Branch("pid",&mPid,"pid/I");
  mTree->Branch("px",&mPx,"px/D");
  mTree->Branch("py",&mPy,"py/D");
  mTree->Branch("pz",&mPz,"pz/D");
  mTree->Branch("vx",&mVx,"vx/D");
  mTree->Branch("vy",&mVy,"vy/D");
  mTree->Branch("vz",&mVz,"vz/D");
  mTree->Branch("theta",&mTheta,"theta/D");
  mTree->Branch("phi",&mPhi,"phi/D");
  mTree->Branch("nhit",&mNHit,"nhit/I");
  mTree->Branch("nhitAegl",&mNHitAegl,"nhitAegl/I");
  mTree->Branch("nhitPhoDet",&mNHitPhoDet,"nhitPhoDet/I");
  mTree->Branch("nelSbKCs",&mNelSbKCs,"nelSbKCs/I");
  mTree->Branch("nelGaAsP",&mNelGaAsP,"nelGaAsP/I");
  mTree->Branch("nelGaAs",&mNelGaAs,"nelGaAs/I");
  mTree->Branch("Lpion",&mLpion,"Lpion/D");
  mTree->Branch("LKaon",&mLKaon,"LKaon/D");
  mTree->Branch("Lproton",&mLproton,"Lproton/D");
  mTree->SetAutoSave(500000);

  return 0;
}

int calLikelihood::Make()
{
  event *aevt = new event(mChainInPut_Events);  /// declear and save info to branchs for event
  hit *ahit = new hit(mChainInPut_Tracks);	/// declear and save info to branchs for track

  long NumOfEvents = (long)mChainInPut_Events->GetEntries();
  mChainInPut_Events->GetEntry(0);
  mChainInPut_Tracks->GetEntry(0);

  // for(int i_event = 0; i_event < 1024; ++i_event) // test event loop
  for(int i_event = 0; i_event < NumOfEvents; ++i_event) // event loop
  { 
    if(i_event%100==0) cout << "processing event:  " << i_event << " ;"<<endl;

    mChainInPut_Events->GetEntry(i_event);  
    mChainInPut_Tracks->GetEntry(i_event);

    const int  pid_gen = aevt->get_pid()->at(0);
    const double px_gen = aevt->get_px()->at(0)/1e3;    //in MeV, convert to GeV
    const double py_gen = aevt->get_py()->at(0)/1e3;    //in MeV, convert to GeV
    const double pz_gen = aevt->get_pz()->at(0)/1e3;    //in MeV, convert to GeV
    const double vx_gen = aevt->get_vx()->at(0);        //in mm
    const double vy_gen = aevt->get_vy()->at(0);        //in mm
    const double vz_gen = aevt->get_vz()->at(0);        //in mm
    const double p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    const double theta_gen=acos(pz_gen/p_gen)*mRICH::DEG;    //in deg
    const double phi_gen=atan2(py_gen,px_gen)*mRICH::DEG;    //in deg

    const int indexSpaceX = utility->get_indexSpaceX(vx_gen);
    const int indexSpaceY = utility->get_indexSpaceY(vy_gen);
    const int indexMomentumP = utility->get_indexMomentumP(px_gen,py_gen,pz_gen);
    const int indexMomentumTheta = utility->get_indexMomentumTheta(px_gen,py_gen,pz_gen);
    const int indexMomentumPhi = utility->get_indexMomentumPhi(px_gen,py_gen);

    int charge = (pid_gen > 0) ? 1 : -1;

    // get pion database
    string key_events_pion = utility->gen_KeyNumOfEvents(charge*211,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_pion = (TH1D*)h_mNumOfEvents[key_events_pion]->Clone();
    string key_photon_pion = utility->gen_KeyMassHypo(charge*211,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_pion = (TH2D*)h_mPhotonDist[key_photon_pion]->Clone();
    int NumofEvents_pion = h_NumofEvents_pion->GetBinContent(1);
    h_database_pion->Sumw2();
    h_database_pion->Scale(1./NumofEvents_pion);

    // get kaon database
    string key_events_kaon = utility->gen_KeyNumOfEvents(charge*321,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_kaon = (TH1D*)h_mNumOfEvents[key_events_kaon]->Clone();
    string key_photon_kaon = utility->gen_KeyMassHypo(charge*321,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_kaon = (TH2D*)h_mPhotonDist[key_photon_kaon]->Clone();
    int NumofEvents_kaon = h_NumofEvents_kaon->GetBinContent(1);
    h_database_kaon->Sumw2();
    h_database_kaon->Scale(1./NumofEvents_kaon);

    // get proton database
    string key_events_proton = utility->gen_KeyNumOfEvents(charge*2212,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_proton = (TH1D*)h_mNumOfEvents[key_events_proton]->Clone();
    string key_photon_proton = utility->gen_KeyMassHypo(charge*2212,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_proton = (TH2D*)h_mPhotonDist[key_photon_proton]->Clone();
    int NumofEvents_proton = h_NumofEvents_proton->GetBinContent(1);
    h_database_proton->Sumw2();
    h_database_proton->Scale(1./NumofEvents_proton);

    // cout << "pid_gen = " << pid_gen << ", charge = " << charge << ", database used: " << endl;
    // cout << key_photon_pion.c_str() << endl;
    // cout << key_photon_kaon.c_str() << endl;
    // cout << key_photon_proton.c_str() << endl;
    // cout << endl;

    // fill photon distribution of unknown PID
    TH2D *h_photonDist_PID = new TH2D("h_photonDist_PID","h_photonDist_PID",mRICH::mNPads,-1.0*mRICH::mHalfWidth,mRICH::mHalfWidth,mRICH::mNPads,-1.0*mRICH::mHalfWidth,mRICH::mHalfWidth);
    int NumOfTracks = ahit->get_hitn()->size();
    for (int i_track = 0; i_track < NumOfTracks; ++i_track) // track loop
    {
      if(isPhoton(ahit,i_track) && !isReflection(ahit,i_track) && isOnPhotonSensor(ahit,i_track))
      {
	double out_x = ahit->get_out_x()->at(i_track);
	double out_y = ahit->get_out_y()->at(i_track);
	h_photonDist_PID->Fill(out_x,out_y);
      }
    }

    mPid = pid_gen;
    mPx  = px_gen;
    mPy  = py_gen;
    mPz  = pz_gen;
    mVx  = vx_gen;
    mVy  = vy_gen;
    mVz  = vz_gen;
    mTheta = theta_gen;
    mPhi = phi_gen;
    mNHit = NumOfTracks;
    mNHitAegl = 0;
    mNHitPhoDet = 0;
    mNelSbKCs = 0;
    mNelGaAsP = 0;
    mNelGaAs = 0;
    mLpion = probability(h_database_pion, h_photonDist_PID);
    mLKaon = probability(h_database_kaon, h_photonDist_PID);
    mLproton = probability(h_database_proton, h_photonDist_PID);

    if(mTree) mTree->Fill();

    delete h_photonDist_PID;
    delete h_NumofEvents_pion;
    delete h_NumofEvents_kaon;
    delete h_NumofEvents_proton;
    delete h_database_pion;
    delete h_database_kaon;
    delete h_database_proton;
  }

  return 0;
}

int calLikelihood::Finish()
{
  cout<<endl;
  cout<<"calLikelihood::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(mFile_OutPut!= NULL){
    mFile_OutPut->cd();
    writeTree();
    mFile_OutPut->Close();
  }
  return 0;
}

int calLikelihood::writeTree()
{
  mTree->Write();
}

bool calLikelihood::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool calLikelihood::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool calLikelihood::isOnAerogel(hit *ahit, int i)
{
  if(ahit->get_out_z()->at(i)>=55.5 && ahit->get_out_z()->at(i)<=85.5) return true;
  //if(ahit->get_out_z()->at(i)>=50.5 && ahit->get_out_z()->at(i)<=70.5) return true;
  else return false;
}

bool calLikelihood::isOnPhotonSensor(hit *ahit, int i)
{
  double out_z = ahit->get_out_z()->at(i);
  if(out_z > 253.6624 && out_z < 255.1626) return true;
  else return false;
}

double calLikelihood::probability(TH2D *h_database, TH2D *h_photonDist_PID)
{
  double prob=1.0;
  const double noise = 2./mRICH::mNPads/mRICH::mNPads; // force every photon from inject particles used in likelihood calculation

  for(unsigned int i_x = 0; i_x < mRICH::mNPads; i_x++)
  {
    for(unsigned int i_y = 0; i_y < mRICH::mNPads; i_y++)
    {
      double k = h_photonDist_PID->GetBinContent(i_x+1,i_y+1); // detected photon number
      double lambda = h_database->GetBinContent(i_x+1,i_y+1); // averaged photon number
      double err_lambda = h_database->GetBinError(i_x+1,i_y+1); // error for specific bin
      if(err_lambda > 0.0) prob *= TMath::PoissonI(k,lambda);
      else prob *= TMath::PoissonI(k,noise);
    }
  }
  return log(prob);
}

////// This is the main function 
int main()
{
  string date = "May10_2018";
  string inputdatabase = Form("/work/eic/xusun/output/database/PDF_database_%s.root",date.c_str());
  string outputfile = Form("/work/eic/xusun/output/likelihood/PID_Likelihood_%s.root",date.c_str());

  calLikelihood *likelihood = new calLikelihood(date,inputdatabase,outputfile);
  
  likelihood->Init();
  likelihood->Make();
  likelihood->Finish();

  return 0;
}
