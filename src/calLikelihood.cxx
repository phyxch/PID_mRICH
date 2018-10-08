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
 mat = new material(); //// initialize the material
 gRandom->SetSeed();
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

  mOutPutFile = Form("/work/eic/xusun/output/likelihood/likelihood_%s_%s.root",mDate.c_str(),mNumOfList.c_str());
  // mOutPutFile = "./out.root"; // batch mode
  cout<<"calLikelihood::Init(), create output file: "<< mOutPutFile.c_str() <<endl;
  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  initChain();
  initHistoMap();
  initHistoQA();
  initGausSmearing();
  initTree();

  return 0;
}

int calLikelihood::initChain()
{
  string inputdir = Form("/work/eic/xusun/output/modular_rich/%s/",mDate.c_str());
  string InPutList = Form("/work/eic/xusun/list/likelihood/%s/mRICH_PID_%s_%s.list",mDate.c_str(),mDate.c_str(),mNumOfList.c_str());
  
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
  for(int i_pid = 0; i_pid < mRICH::mNumOfParType; ++i_pid)
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
  h_mPhotonEnergy_Wavelength = new TH2D("h_mPhotonEnergy_Wavelength","h_mPhotonEnergy_Wavelength",100,0.0,10.0,1000,0.0,1000.0);
  hnHitAeglPerEvtvsMom = new TH2D ("hnHitAeglPerEvtvsMom","hnHitAeglPerEvtvsMom",200,0.,20.0,400,0.,400.);
  hnPhotonElvsnHits_SbKCs = new TH2D ("hnPhotonElvsnHits_SbKCs","hnPhotonElvsnHits_SbKCs",100,0.,100.,50,0.,50.);
  hnPhotonElvsnHits_GaAsP = new TH2D ("hnPhotonElvsnHits_GaAsP","hnPhotonElvsnHits_GaAsP",100,0.,100.,50,0.,50.);
  hnPhotonElvsnHits_GaAs = new TH2D ("hnPhotonElvsnHits_GaAs","hnPhotonElvsnHits_GaAs",100,0.,100.,50,0.,50.);

  h_mXGausSmearing = new TH2D("h_mXGausSmearing","h_mXGausSmearing",121,-60.5,60.5,121,-60.5,60.5);
  h_mYGausSmearing = new TH2D("h_mYGausSmearing","h_mYGausSmearing",121,-60.5,60.5,121,-60.5,60.5);

  return 0;
}

int calLikelihood::initGausSmearing()
{
  f_mGaus = new TF1("f_mGaus","gaus",-20.0,20.0);
  f_mGaus->SetParameter(0,1.0);
  f_mGaus->SetParameter(1,0.0);
  f_mGaus->SetParameter(2,1.0);
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
  mTree->Branch("Lelectron",&mLelectron,"Lelectron/D");
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

    // int charge = (pid_gen > 0) ? 1 : -1;
    int charge = -1; // temperary for electron implementation

    // get electron database
    string key_events_electron = utility->gen_KeyNumOfEvents(11,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_electron = (TH1D*)h_mNumOfEvents[key_events_electron]->Clone();
    string key_photon_electron = utility->gen_KeyMassHypo(11,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_electron = (TH2D*)h_mPhotonDist[key_photon_electron]->Clone();
    int NumofEvents_electron = h_NumofEvents_electron->GetBinContent(1);
    h_database_electron->Sumw2();
    if(NumofEvents_electron > 0) h_database_electron->Scale(1./NumofEvents_electron);


    // get pion database
    string key_events_pion = utility->gen_KeyNumOfEvents(charge*211,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_pion = (TH1D*)h_mNumOfEvents[key_events_pion]->Clone();
    string key_photon_pion = utility->gen_KeyMassHypo(charge*211,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_pion = (TH2D*)h_mPhotonDist[key_photon_pion]->Clone();
    int NumofEvents_pion = h_NumofEvents_pion->GetBinContent(1);
    h_database_pion->Sumw2();
    if(NumofEvents_pion > 0) h_database_pion->Scale(1./NumofEvents_pion);

    // get kaon database
    string key_events_kaon = utility->gen_KeyNumOfEvents(charge*321,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_kaon = (TH1D*)h_mNumOfEvents[key_events_kaon]->Clone();
    string key_photon_kaon = utility->gen_KeyMassHypo(charge*321,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_kaon = (TH2D*)h_mPhotonDist[key_photon_kaon]->Clone();
    int NumofEvents_kaon = h_NumofEvents_kaon->GetBinContent(1);
    h_database_kaon->Sumw2();
    if(NumofEvents_kaon > 0) h_database_kaon->Scale(1./NumofEvents_kaon);

    // get proton database
    string key_events_proton = utility->gen_KeyNumOfEvents(charge*2212,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH1D *h_NumofEvents_proton = (TH1D*)h_mNumOfEvents[key_events_proton]->Clone();
    string key_photon_proton = utility->gen_KeyMassHypo(charge*2212,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    TH2D *h_database_proton = (TH2D*)h_mPhotonDist[key_photon_proton]->Clone();
    int NumofEvents_proton = h_NumofEvents_proton->GetBinContent(1);
    h_database_proton->Sumw2();
    if(NumofEvents_proton > 0) h_database_proton->Scale(1./NumofEvents_proton);

    // cout << "pid_gen = " << pid_gen << ", charge = " << charge << ", database used: " << endl;
    // cout << key_photon_pion.c_str() << endl;
    // cout << key_photon_kaon.c_str() << endl;
    // cout << key_photon_proton.c_str() << endl;
    // cout << endl;

    // fill photon distribution of unknown PID
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

	// double QE_GaAsP = mat->extrapQE_GaAsP(wavelength); // quantum efficiency of photon sensor => need to be updated
	// if(QE_GaAsP > gRandom->Uniform(0.0,1.0))
	double QE = mat->extrapQE(wavelength); // get quantum efficiency for photon sensor => need to be updated
	// cout << "wavelength = " << wavelength << ", QE_GaAsP = " << QE_GaAsP << ", QE = " << QE << endl;
	if( QE > gRandom->Uniform(0.0,1.0) )
	{
	  NumOfElGaAsP++; 
	  double out_x_input = ahit->get_out_x()->at(i_track);
	  double out_y_input = ahit->get_out_y()->at(i_track);
	  double delta_x = GausSmearing(f_mGaus);
	  double delta_y = GausSmearing(f_mGaus);

	  double out_x = out_x_input+delta_x;
	  double out_y = out_y_input+delta_y;
	  if( isInSensorPlane(out_x,out_y) )
	  {
	    h_photonDist_PID->Fill(out_x,out_y);
	    // cout << "out_x_input = " << out_x_input << ", out_x = " << out_x << endl;
	    // cout << "out_y_input = " << out_y_input << ", out_y = " << out_y << endl;
	    // cout << endl;
	  }
	  h_mXGausSmearing->Fill(out_x_input,out_x);
	  h_mYGausSmearing->Fill(out_y_input,out_y);
	}

	double QE_GaAs = mat->extrapQE_GaAs(wavelength);
	if(QE_GaAs > gRandom->Uniform(0.0,1.0))
	{
	  NumOfElGaAs++;
	}
      }
    }

    // add noise with 2 electrons
    for(int i_electron = 0; i_electron < 2; i_electron++) 
    {
      double out_x = gRandom->Uniform(2.5,mRICH::mHalfWidth-2.0);
      double out_y = gRandom->Uniform(2.5,mRICH::mHalfWidth-2.0);
      double sign_x = 1.0;
      double sign_y = 1.0;
      if(gRandom->Rndm() > 0.5) sign_x = -1.0;
      if(gRandom->Rndm() > 0.5) sign_y = -1.0;
      h_photonDist_PID->Fill(out_x*sign_x,out_y*sign_y);
      // cout << "out_x = " << out_x << ", sign_x = " << sign_x << endl;
      // cout << "out_y = " << out_y << ", sign_y = " << sign_y << endl;
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
    mLelectron = probability(h_database_electron, h_photonDist_PID);

    if(mTree) mTree->Fill();

    delete h_photonDist_PID;
    delete h_NumofEvents_pion;
    delete h_NumofEvents_kaon;
    delete h_NumofEvents_proton;
    delete h_database_pion;
    delete h_database_kaon;
    delete h_database_proton;
    delete h_database_electron;
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
    writeQA();
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
  h_mXGausSmearing->Write();
  h_mYGausSmearing->Write();

  return 0;
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
  // if(ahit->get_out_z()->at(i)>=63.5874 && ahit->get_out_z()->at(i)<=96.5876) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 1) return true;
  else return false;
}

bool calLikelihood::isOnPhotonSensor(hit *ahit, int i)
{
  // double out_z = ahit->get_out_z()->at(i);
  // if(out_z > 253.6624+3.0 && out_z < 255.1626+3.0) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 2) return true;
  else return false;
}

double calLikelihood::probability(TH2D *h_database, TH2D *h_photonDist_PID)
{
  double prob=1.0;
  const double noise = 2.0/mRICH::mNumOfPixels/mRICH::mNumOfPixels; // force every photon from inject particles used in likelihood calculation

  for(unsigned int i_x = 0; i_x < mRICH::mNumOfPixels; i_x++)
  {
    for(unsigned int i_y = 0; i_y < mRICH::mNumOfPixels; i_y++)
    {
      double k = h_photonDist_PID->GetBinContent(i_x+1,i_y+1); // detected photon number
      double err_k = h_photonDist_PID->GetBinError(i_x+1,i_y+1); // error for specific bin
      double lambda = h_database->GetBinContent(i_x+1,i_y+1); // averaged photon number
      double err_lambda = h_database->GetBinError(i_x+1,i_y+1); // error for specific bin
      if(err_k > 0)
      {
	if(err_lambda > 0.0) prob *= TMath::PoissonI(k,lambda);
	else prob *= TMath::PoissonI(k,noise);
      }
    }
  }
  return log(prob);
}

double calLikelihood::GausSmearing(TF1 *f_gaus)
{
  double delta_pos = f_gaus->GetRandom();
  return delta_pos;
}

bool calLikelihood::isInSensorPlane(double out_x, double out_y)
{
  if( !(TMath::Abs(out_x) >= 2.5 && TMath::Abs(out_x) <= mRICH::mHalfWidth-2.0) ) return false;
  if( !(TMath::Abs(out_y) >= 2.5 && TMath::Abs(out_y) <= mRICH::mHalfWidth-2.0) ) return false;
  return true;
}

////// This is the main function 
int main(int argc, char **argv)
{
  if(argc!=2) return 0;

  const char *input = argv[1];
  string numoflist(input);
  
  string date = "Oct08_2018";
  string inputdatabase = Form("/work/eic/xusun/output/database/database_%s.root",date.c_str());

  cout << "numoflist = " << numoflist.c_str() << endl;
  calLikelihood *likelihood = new calLikelihood(numoflist,date,inputdatabase);
  
  likelihood->Init();
  likelihood->Make();
  likelihood->Finish();

  cout << "This is the end of calLikelihood!!!" << endl;

  return 0;
}
