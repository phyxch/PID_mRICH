#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include <utility>
#include "string.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TString.h"
#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/genMassHypo.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;
using namespace TMath;


genMassHypo::genMassHypo(string numoflist, string date)
{
 cout<<endl;
 cout<<"genMassHypo::genMassHypo() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mNumOfList = numoflist;
 mDate = date;
 utility = new Utility(); // initialize utility class
 mat = new material(); //// initialize the material
 gRandom->SetSeed();
}

genMassHypo::~genMassHypo()
{
 cout<<"genMassHypo::~genMassHypo() ----- Release memory ! ------"<<endl;
 delete mat;
 delete utility;
 delete File_mOutPut;
}

int genMassHypo::Init()
{
  cout<<"genMassHypo::Init() ----- Initialization ! ------"<<endl;

  // mOutPutFile = Form("/work/eic/xusun/output/database/database_%s_%s.root",mDate.c_str(),mNumOfList.c_str());
  mOutPutFile = "./out.root"; // batch mode
  cout<<"genMassHypo::Init(), create output file: "<< mOutPutFile.c_str() <<endl;
  File_mOutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  initChain();
  initHistoMap();
  initGausSmearing();
  initHistoQA();

  return 0;
}

int genMassHypo::initChain()
{
  string inputdir = Form("/work/eic/xusun/output/modular_rich/%s/",mDate.c_str());
  // string inputdir = "/work/eic/xusun/output/modular_rich/BeamTest/PlusMinus/";
  string InPutList = Form("/work/eic/xusun/list/database/%s/mRICH_PDF_%s_%s.list",mDate.c_str(),mDate.c_str(),mNumOfList.c_str());
  
  mChainInPut_Events = new TChain("generated");
  mChainInPut_Tracks = new TChain("eic_rich");

  if (!InPutList.empty())   // if input file is ok
  {
    cout << "Open input database file list: " << InPutList.c_str() << endl;
    ifstream in(InPutList.c_str());  // input stream
    if(in)
    {
      cout << "input database file list is ok" << endl;
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
      cout << "WARNING: input database file input is problemtic" << endl;
    }
  }

  long NumOfEvents = (long)mChainInPut_Events->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;

  return 0;
}

int genMassHypo::initHistoMap()
{
  cout<<"genMassHypo::initHistoMap(), initialize database histograms ;"<<endl;

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
	      // cout << "genMassHypo::init(), initialize histogram: " << key_events.c_str() << endl;
	      h_mNumOfEvents[key_events] = new TH1D(key_events.c_str(),key_events.c_str(),1,-0.5,0.5);

	      string key_photon = utility->gen_KeyMassHypo(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "genMassHypo::init(), initialize histogram: " << key_photon.c_str() << endl;
	      h_mPhotonDist[key_photon] = new TH2D(key_photon.c_str(),key_photon.c_str(),mRICH::mNumOfPixels,mRICH::mPixels,mRICH::mNumOfPixels,mRICH::mPixels);

	      string key_photon_generated = utility->gen_KeyMassHypoGenerated(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "genMassHypo::init(), initialize histogram: " << key_photon_generated.c_str() << endl;
	      h_mPhotonGenerated[key_photon_generated] = new TH2D(key_photon_generated.c_str(),key_photon_generated.c_str(),mRICH::mNumOfPixels,mRICH::mPixels,mRICH::mNumOfPixels,mRICH::mPixels);
	    }
	  }
	}
      }
    }
  }

  return 0;
}

int genMassHypo::initGausSmearing()
{
  f_mGaus = new TF1("f_mGaus","gaus",-20.0,20.0);
  f_mGaus->SetParameter(0,1.0);
  f_mGaus->SetParameter(1,0.0);
  f_mGaus->SetParameter(2,1.0);

  return 0;
}

int genMassHypo::initHistoQA()
{
  h_mXGausSmearing = new TH2D("h_mXGausSmearing","h_mXGausSmearing",121,-60.5,60.5,121,-60.5,60.5);
  h_mYGausSmearing = new TH2D("h_mYGausSmearing","h_mYGausSmearing",121,-60.5,60.5,121,-60.5,60.5);

  return 0;
}

int genMassHypo::Make()
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

    const double momentum = TMath::Sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
    const double theta = TMath::ACos(pz_gen/momentum)*mRICH::DEG;    //in deg
    const double phi = TMath::ATan2(py_gen,px_gen)*mRICH::DEG;    //in deg            

    const int indexSpaceX = utility->get_indexSpaceX(vx_gen);
    const int indexSpaceY = utility->get_indexSpaceY(vy_gen);
    const int indexMomentumP = utility->get_indexMomentumP(px_gen,py_gen,pz_gen);
    const int indexMomentumTheta = utility->get_indexMomentumTheta(px_gen,py_gen,pz_gen);
    const int indexMomentumPhi = utility->get_indexMomentumPhi(px_gen,py_gen);

    if(indexSpaceX < 0 || indexSpaceY < 0 || indexMomentumP < 0 || indexMomentumTheta < 0 || indexMomentumPhi < 0) continue;

    // cout << "indexMomentumP = " << indexMomentumP << ", indexMomentumTheta = " << indexMomentumTheta << ", indexMomentumPhi = " << indexMomentumPhi << endl;

    string key_events = utility->gen_KeyNumOfEvents(pid_gen,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    // cout << "fill histogram: " << key_events.c_str() << endl;
    h_mNumOfEvents[key_events]->Fill(0);

    string key_photon = utility->gen_KeyMassHypo(pid_gen,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    // cout << "fill histogram: " << key_photon.c_str() << endl;

    string key_photon_generated = utility->gen_KeyMassHypoGenerated(pid_gen,indexSpaceX,indexSpaceY,indexMomentumP,indexMomentumTheta,indexMomentumPhi);
    // cout << "fill histogram: " << key_photon_generated.c_str() << endl;

    int NumOfTracks = ahit->get_hitn()->size();
    for (int i_track = 0; i_track < NumOfTracks; ++i_track) // track loop
    {
      if(isPhoton(ahit,i_track) && !isReflection(ahit,i_track) && isOnPhotonSensor(ahit,i_track))
      {
	double out_x_generated = ahit->get_out_x()->at(i_track);
	double out_y_generated = ahit->get_out_y()->at(i_track);
	h_mPhotonGenerated[key_photon_generated]->Fill(out_x_generated,out_y_generated);

	double photonE = ahit->get_trackE()->at(i_track);   /// in MeV (GEANT4 default)
	double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
	double QE_GaAsP = mat->extrapQE_GaAsP(wavelength); // get quantum efficiency for photon sensor => need to be updated

	// if( QE_GaAsP > gRandom->Uniform(0.0,1.0) )
	{
	  double out_x_input = ahit->get_out_x()->at(i_track);
	  double out_y_input = ahit->get_out_y()->at(i_track);
	  double delta_x = GausSmearing(f_mGaus);
	  double delta_y = GausSmearing(f_mGaus);

	  double out_x = out_x_input+delta_x;
	  double out_y = out_y_input+delta_y;
	  if( isInSensorPlane(out_x,out_y) )
	  {
	    h_mPhotonDist[key_photon]->Fill(out_x,out_y);
	    // cout << "out_x_input = " << out_x_input << ", out_x = " << out_x << endl;
	    // cout << "out_y_input = " << out_y_input << ", out_y = " << out_y << endl;
	    // cout << endl;
	  }
	  h_mXGausSmearing->Fill(out_x_input,out_x);
	  h_mYGausSmearing->Fill(out_y_input,out_y);
	}
      }
    }
  }

  return 0;
}

int genMassHypo::Finish()
{
  cout<<endl;
  cout<<"genMassHypo::Finish() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(File_mOutPut != NULL){
    File_mOutPut->cd();
    writeHistoMap();
    writeHistoQA();
    File_mOutPut->Close();
  }
  return 0;
}

int genMassHypo::writeHistoMap()
{
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
	      // cout << "genMassHypo::end(), write histogram: " << key_events.c_str() << endl;
	      h_mNumOfEvents[key_events]->Write();

	      string key_photon = utility->gen_KeyMassHypo(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "genMassHypo::end(), write histogram: " << key_photon.c_str() << endl;
	      h_mPhotonDist[key_photon]->Write();

	      string key_photon_generated = utility->gen_KeyMassHypoGenerated(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_mom,i_theta,i_phi);
	      // cout << "genMassHypo::end(), write histogram: " << key_photon_generated.c_str() << endl;
	      h_mPhotonGenerated[key_photon_generated]->Write();
	    }
	  }
	}
      }
    }
  }

  return 0;
}

int genMassHypo::writeHistoQA()
{
  h_mXGausSmearing->Write();
  h_mYGausSmearing->Write();
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
  // if(ahit->get_out_z()->at(i)>=63.5874 && ahit->get_out_z()->at(i)<=96.5876) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 1) return true;
  else return false;
}

bool genMassHypo::isOnPhotonSensor(hit *ahit, int i)
{
  // double out_z = ahit->get_out_z()->at(i);
  // if(out_z > 253.6624+3.0 && out_z < 255.1626+3.0) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 2) return true;
  else return false;
}

double genMassHypo::GausSmearing(TF1 *f_gaus)
{
  double delta_pos = f_gaus->GetRandom();
  return delta_pos;
}

bool genMassHypo::isInSensorPlane(double out_x, double out_y)
{
  if( !(TMath::Abs(out_x) > 2.5 && TMath::Abs(out_x) < mRICH::mHalfWidth-2.0) ) return false;
  if( !(TMath::Abs(out_y) > 2.5 && TMath::Abs(out_y) < mRICH::mHalfWidth-2.0) ) return false;
  return true;
}

////// This is the main function 
int main(int argc, char **argv)
{
  if(argc!=2) return 0;

  const char *input = argv[1];
  string numoflist(input);
  
  string date = "Oct08_2018";
  
  cout << "numoflist = " << numoflist.c_str() << endl;
  genMassHypo *genMassHypotheses = new genMassHypo(numoflist,date);
  
  genMassHypotheses->Init();
  genMassHypotheses->Make();
  genMassHypotheses->Finish();

  cout << "This is the end of genMassHypo!!!" << endl;
  
  return 0;
}
