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
#include <TRandom3.h>
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/calLikelihood.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;

calLikelihood::calLikelihood(string numoflist, string date, string inputdatabase)
{
 cout<<endl;
 cout<<"calLikelihood::calLikelihood() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mInPutDataBase = inputdatabase;
 mNumOfList = numoflist;
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

  // mOutPutFile = Form("/work/eic/xusun/output/likelihood/likelihood_%s_%s.root",mDate.c_str(),mNumOfList.c_str());
  mOutPutFile = "./out.root"; // batch mode
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
  string InPutList = Form("/work/eic/xusun/list/likelihood/mRICH_PID_%s_%s.list",mDate.c_str(),mNumOfList.c_str());
  
  mChainInPut_Events = new TChain("generated");
  mChainInPut_Tracks = new TChain("eic_rich");

  if (!InPutList.empty())   // if input file is ok
  {
    cout << "Open input likelihood file list " << endl;
    ifstream in(InPutList.c_str());  // input stream
    if(in)
    {
      cout << "input likelihood file list is ok" << endl;
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
      cout << "WARNING: input likelihood file input is problemtic" << endl;
    }
  }

  long NumOfEvents = (long)mChainInPut_Events->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;

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

int calLikelihood::initHistoQA()
{
  h_mPhotonEnergy_Wavelength = new TH2D("h_mPhotonEnergy_Wavelength","h_mPhotonEnergy_Wavelength",100,0.0,5.0,100,250.0,750.0);
  hnHitAeglPerEvtvsMom = new TH2D ("hnHitAeglPerEvtvsMom","hnHitAeglPerEvtvsMom",100,0.,10.,200,0.,200.);
  hnPhotonElvsnHits_SbKCs = new TH2D ("hnPhotonElvsnHits_SbKCs","hnPhotonElvsnHits_SbKCs",50,0.,50.,25,0.,25.);
  hnPhotonElvsnHits_GaAsP = new TH2D ("hnPhotonElvsnHits_GaAsP","hnPhotonElvsnHits_GaAsP",50,0.,50.,25,0.,25.);
  hnPhotonElvsnHits_GaAs = new TH2D ("hnPhotonElvsnHits_GaAs","hnPhotonElvsnHits_GaAs",50,0.,50.,25,0.,25.);
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
  mTree->Branch("nelSbKCs",&mNElSbKCs,"nelSbKCs/I");
  mTree->Branch("nelGaAsP",&mNElGaAsP,"nelGaAsP/I");
  mTree->Branch("nelGaAs",&mNElGaAs,"nelGaAs/I");
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
    if(i_event%1000==0) cout << "processing event:  " << i_event << " ;"<<endl;

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

    if(indexSpaceX < 0 || indexSpaceY < 0 || indexMomentumP < 0 || indexMomentumTheta < 0 || indexMomentumPhi < 0) continue;

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
    // TH2D *h_photonDist_PID = new TH2D("h_photonDist_PID","h_photonDist_PID",mRICH::mNumOfPixels,-1.0*mRICH::mHalfWidth,mRICH::mHalfWidth,mRICH::mNumOfPixels,-1.0*mRICH::mHalfWidth,mRICH::mHalfWidth);
    TH2D *h_photonDist_PID = new TH2D("h_photonDist_PID","h_photonDist_PID",mRICH::mNumOfPixels,mRICH::mPixels,mRICH::mNumOfPixels,mRICH::mPixels);
    int NumOfHitAegl = 0;
    int NumOfHitPhoDet = 0;
    int NumOfElSbKCs = 0;
    int NumOfElGaAsP = 0;
    int NumOfElGaAs = 0;
    int NumOfTracks = ahit->get_hitn()->size();
    for (int i_track = 0; i_track < NumOfTracks; ++i_track) // track loop
    {
      double photonE = ahit->get_trackE()->at(i_track);   // in MeV (GEANT4 default)
      double wavelength = 1240.0/(photonE*1.0e6);  // MeV->eV,wavelength in "nm"

      if(isPhoton(ahit,i_track) && !isReflection(ahit,i_track) && isOnAerogel(ahit,i_track))
      {
	NumOfHitAegl++;
	h_mPhotonEnergy_Wavelength->Fill(photonE*1.0e6,wavelength);
      }

      if(isPhoton(ahit,i_track) && !isReflection(ahit,i_track) && isOnPhotonSensor(ahit,i_track))
      {
	NumOfHitPhoDet++;

	double QE_SbKCs = mat->extrapQE_SbKCs(wavelength);
	if(QE_SbKCs > gRandom->Uniform(0.0,1.0))
	{
	  NumOfElSbKCs++;
	}

	double QE_GaAsP = mat->extrapQE_GaAsP(wavelength); // quantum efficiency of photon sensor => need to be updated
	if(QE_GaAsP > gRandom->Uniform(0.0,1.0))
	{
	  NumOfElGaAsP++; 
	  double out_x = ahit->get_out_x()->at(i_track);
	  double out_y = ahit->get_out_y()->at(i_track);
	  h_photonDist_PID->Fill(out_x,out_y);
	}

	double QE_GaAs = mat->extrapQE_GaAs(wavelength);
	if(QE_GaAs > gRandom->Uniform(0.0,1.0))
	{
	  NumOfElGaAs++;
	}
      }
    }

    for(int i_electron = 0; i_electron < 2; i_electron++) 
    {
      double out_x = gRandom->Uniform(-1.0,1.0)*mRICH::mHalfWidth;
      double out_y = gRandom->Uniform(-1.0,1.0)*mRICH::mHalfWidth;
      h_photonDist_PID->Fill(out_x,out_y);
    }
  

    hnHitAeglPerEvtvsMom->Fill(p_gen,NumOfHitAegl);
    hnPhotonElvsnHits_SbKCs->Fill(NumOfHitPhoDet,NumOfElSbKCs);
    hnPhotonElvsnHits_GaAsP->Fill(NumOfHitPhoDet,NumOfElGaAsP);
    hnPhotonElvsnHits_GaAs->Fill(NumOfHitPhoDet,NumOfElGaAs);

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
    mNHitAegl = NumOfHitAegl;
    mNHitPhoDet = NumOfHitPhoDet;
    mNElSbKCs = NumOfElSbKCs;
    mNElGaAsP = NumOfElGaAsP;
    mNElGaAs = NumOfElGaAs;
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

  return 0;
}

int calLikelihood::writeQA()
{
  h_mPhotonEnergy_Wavelength->Write();
  hnHitAeglPerEvtvsMom->Write();
  hnPhotonElvsnHits_SbKCs->Write();
  hnPhotonElvsnHits_GaAsP->Write();
  hnPhotonElvsnHits_GaAs->Write();
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
  const double noise = 2./mRICH::mNumOfPixels/mRICH::mNumOfPixels; // force every photon from inject particles used in likelihood calculation

  for(unsigned int i_x = 0; i_x < mRICH::mNumOfPixels; i_x++)
  {
    for(unsigned int i_y = 0; i_y < mRICH::mNumOfPixels; i_y++)
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
int main(int argc, char **argv)
{
  if(argc!=2) return 0;

  const char *input = argv[1];
  string numoflist(input);
  
  string date = "May23_2018";
  string inputdatabase = Form("/work/eic/xusun/output/database/database_%s.root",date.c_str());

  cout << "numoflist = " << numoflist.c_str() << endl;
  calLikelihood *likelihood = new calLikelihood(numoflist,date,inputdatabase);
  
  likelihood->Init();
  likelihood->Make();
  likelihood->Finish();

  cout << "This is the end of calLikelihood!!!" << endl;

  return 0;
}
