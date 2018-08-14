#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include <utility>
#include "string.h"

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"

#include "../include/event.h"
#include "../include/hit.h"
#include "../include/material.h"
#include "../include/RingFinder.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;
using namespace TMath;


RingFinder::RingFinder(string numoflist, string date)
{
 cout<<endl;
 cout<<"RingFinder::RingFinder() ----- Constructor ! ------"<<endl;
 cout<<endl;
 mNumOfList = numoflist;
 mDate = date;
 utility = new Utility(); // initialize utility class
 mat = new material(); //// initialize the material
 gRandom->SetSeed();
}

RingFinder::~RingFinder()
{
 cout<<"RingFinder::~RingFinder() ----- Release memory ! ------"<<endl;
 delete mat;
 delete utility;
 delete File_mOutPut;
}

int RingFinder::Init()
{
  cout<<"RingFinder::Init() ----- Initialization ! ------"<<endl;

  mOutPutFile = Form("/work/eic/xusun/output/database/database_%s_%s.root",mDate.c_str(),mNumOfList.c_str());
  // mOutPutFile = "./out.root"; // batch mode
  cout<<"RingFinder::Init(), create output file: "<< mOutPutFile.c_str() <<endl;
  File_mOutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  initChain();
  initHistoMap();
  initGausSmearing();
  initHistoQA();

  return 0;
}

int RingFinder::initChain()
{
  // string inputdir = Form("/work/eic/xusun/output/modular_rich/%s/",mDate.c_str());
  string inputdir = "/work/eic/xusun/output/modular_rich/BeamTest/PlusMinus/";
  string InPutList = Form("/work/eic/xusun/list/database/mRICH_PDF_%s_%s.list",mDate.c_str(),mNumOfList.c_str());
  
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

int RingFinder::initHistoMap()
{
  cout<<"RingFinder::initHistoMap(), initialize Ring Finder histogram;"<<endl;

  h_mPhotonGenerated = new TH2D("h_mPhotonGenerated","h_mPhotonGenerated",mRICH::mNumOfPixels,mRICH::mPixels,mRICH::mNumOfPixels,mRICH::mPixels);

  h_mNumOfPhotons = new TH1D("h_mNumOfPhotons","h_mNumOfPhotons",1,-0.5,0.5);
  h_mPhotonDist = new TH2D("h_mPhotonDist","h_mPhotonDist",mRICH::mNumOfPixels,mRICH::mPixels,mRICH::mNumOfPixels,mRICH::mPixels);

  clearHistoMap();

  return 0;
}

int RingFinder::initGausSmearing()
{
  f_mGaus = new TF1("f_mGaus","gaus",-20.0,20.0);
  f_mGaus->SetParameter(0,1.0);
  f_mGaus->SetParameter(1,0.0);
  f_mGaus->SetParameter(2,1.0);

  return 0;
}

int RingFinder::initHistoQA()
{
  h_mXGausSmearing = new TH2D("h_mXGausSmearing","h_mXGausSmearing",121,-60.5,60.5,121,-60.5,60.5);
  h_mYGausSmearing = new TH2D("h_mYGausSmearing","h_mYGausSmearing",121,-60.5,60.5,121,-60.5,60.5);

  return 0;
}

//--------------------------------------------------------

int RingFinder::Make()
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

    // const int indexSpaceX = utility->get_indexSpaceX(vx_gen);
    // const int indexSpaceY = utility->get_indexSpaceY(vy_gen);
    // const int indexMomentumP = utility->get_indexMomentumP(px_gen,py_gen,pz_gen);
    // const int indexMomentumTheta = utility->get_indexMomentumTheta(px_gen,py_gen,pz_gen);
    // const int indexMomentumPhi = utility->get_indexMomentumPhi(px_gen,py_gen);
    //
    // if(indexSpaceX < 0 || indexSpaceY < 0 || indexMomentumP < 0 || indexMomentumTheta < 0 || indexMomentumPhi < 0) continue;


    int NumOfTracks = ahit->get_hitn()->size();
    for (int i_track = 0; i_track < NumOfTracks; ++i_track) // track loop
    {
      if(isPhoton(ahit,i_track) && !isReflection(ahit,i_track) && isOnPhotonSensor(ahit,i_track))
      {
	double out_x_generated = ahit->get_out_x()->at(i_track);
	double out_y_generated = ahit->get_out_y()->at(i_track);
	h_mPhotonGenerated->Fill(out_x_generated,out_y_generated);

	double photonE = ahit->get_trackE()->at(i_track);   /// in MeV (GEANT4 default)
	double wavelength = 1240./(photonE*1.e6);  /// MeV->eV,wavelength in "nm"
	double QE_GaAsP = mat->extrapQE_GaAsP(wavelength); // get quantum efficiency for photon sensor => need to be updated

	if( QE_GaAsP > gRandom->Uniform(0.0,1.0) )
	{
	  double out_x_input = ahit->get_out_x()->at(i_track);
	  double out_y_input = ahit->get_out_y()->at(i_track);
	  double delta_x = GausSmearing(f_mGaus);
	  double delta_y = GausSmearing(f_mGaus);

	  double out_x = out_x_input+delta_x;
	  double out_y = out_y_input+delta_y;
	  if( isInSensorPlane(out_x,out_y) )
	  {
	    h_mNumOfPhotons->Fill(0);
	    h_mPhotonDist->Fill(out_x,out_y);
	    int binX = h_mPhotonDist->GetXais()->FindBin(out_x);
	    int binY = h_mPhotonDist->GetYais()->FindBin(out_y);
	    mXPixelMap.push_back(binX);
	    mYPixelMap.push_back(binY);
	    // cout << "out_x_input = " << out_x_input << ", out_x = " << out_x << endl;
	    // cout << "out_y_input = " << out_y_input << ", out_y = " << out_y << endl;
	    // cout << endl;
	  }
	  h_mXGausSmearing->Fill(out_x_input,out_x);
	  h_mYGausSmearing->Fill(out_y_input,out_y);
	}
      }
    }
    HoughTransform(h_mNumOfPhotons,h_mPhotonDist, mXPixelMap, mYPixelMap);
    clearHistoMap();
  }

  return 0;
}

int RingFinder::clearHistoMap()
{
  h_mPhotonGenerated->Reset();

  h_mNumOfPhotons->Reset();
  h_mPhotonDist->Reset();
  mXPixelMap.clear();
  mYPixelMap.clear();

  return 0;
}

int RingFinder::HoughTransform(TH1D *h_NumOfPhotons, TH2D *h_PhotonDist, intVec xPixel, intVec yPixel);
{
  // int NumOfPhotons = h_NumOfPhotons->GetBinContent(1);
  int NumOfPhotons = h_NumOfPhotons->GetEntries();
  for(int i_hit_1st = 0; i_hit_1st < NumOfPhotons-2; ++i_hit_1st)
  {
    hitPosition firstHit;
    firstHit.x = h_PhotonDist->GetXais()->GetBinContent(xPixel[i_hit_1st]);
    firstHit.y = h_PhotonDist->GetYais()->GetBinContent(yPixel[i_hit_1st]);
    for(int i_hit_2nd = 1; i_hit_2nd < NumOfPhotons-1; ++i_hit_2nd)
    {
      hitPosition secondHit;
      secondHit.x = h_PhotonDist->GetXais()->GetBinContent(xPixel[i_hit_2nd]);
      secondHit.y = h_PhotonDist->GetYais()->GetBinContent(yPixel[i_hit_2nd]);
      for(int i_hit_3rd = 2; i_hit_3rd < NumOfPhotons; ++i_hit_3rd)
      {
	hitPosition thirdHit;
	thirdHit.x = h_PhotonDist->GetXais()->GetBinContent(xPixel[i_hit_3rd]);
	thirdHit.y = h_PhotonDist->GetYais()->GetBinContent(yPixel[i_hit_3rd]);
	double x_Cherenkov = -999.9
	double y_Cherenkov = -999.9
	double r_Cherenkov = -999.9

	bool ringStatus = findRing(firstHit,secondHit,thirdHit, x_Cherenkov, y_Cherenkov, r_Cherenkov);
	if(ringStatus) 
	{
	}
      }
    }
  }
}

bool RingFinder::findRing(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit, double &x_Cherenkov, double &y_Cherenkov, double &r_Cherenkov)
{
  // check if 3 hit points are at same position or collinear
  if(isSamePosition(firstHit,secondHit,thirdHit) || isCollinear(firstHit,secondHit,thirdHit) ) return false;

  // calculate middle point between 1st $ 2nd and 1st & 3rd
  hitPosition midPoint12, midPoint13;
  midPoint12.x = (firstHit.x+secondHit.x)/2.0;
  midPoint12.y = (firstHit.y+secondHit.y)/2.0;
  double slope12 = -(firstHit.x-secondHit.x)/(firstHit.y-secondHit.y);

  midPoint13.x = (firstHit.x+thirdHit.x)/2.0;
  midPoint13.y = (firstHit.y+thirdHit.y)/2.0;
  double slope13 = -(firstHit.x-thirdHit.x)/(firstHit.y-thirdHit.y);

  // y_Cherenkov - midPoint12.y = slope12*(x_Cherenkov - midPoint12.x)
  x_Cherenkov = ((slope12*midPoint12.x-midPoint12.y)-(slope13*midPoint13.x-midPoint13.y))/(slope12-slope13);
  y_Cherenkov = midPoint12.y+slope12*(x_Cherenkov-midPoint12.x);
  // y_Cherenkov = (slope13*(midPoint12.y-slope12*midPoint12.x)-slope12*(midPoint13.y-slope13*midPoint13.x))/(slope13-slope12);
  r_Cherenkov = TMath::Sqrt((x_Cherenkov-firstHit.x)*(x_Cherenkov-firstHit.x)+(y_Cherenkov-firstHit.y)*(y_Cherenkov-firstHit.y));

  return true;
}

bool RingFinder::isSamePosition(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit)
{
  if(TMath::Abs(firstHit.x-secondHit.x) < 1e-5 && TMath::Abs(firstHit.y-secondHit.y) < 1e-5) return true;
  if(TMath::Abs(firstHit.x-thirdHit.x)  < 1e-5 && TMath::Abs(firstHit.y-thirdHit.y)  < 1e-5) return true;
  if(TMath::Abs(secondHit.x-thirdHit.x) < 1e-5 && TMath::Abs(secondHit.y-thirdHit.y) < 1e-5) return true;

  return false;
}

bool RingFinder::isCollinear(hitPosition firstHit, hitPosition secondHit, hitPosition thirdHit)
{
  if(TMath::Abs(firstHit.x-secondHit.x) < 1e-5 && TMath::Abs(firstHit.x-thirdHit.x) < 1e-5) return true;
  if(TMath::Abs(firstHit.y-secondHit.y) < 1e-5 && TMath::Abs(firstHit.y-thirdHit.y) < 1e-5) return true;

  double slope12 = (firstHit.y-secondHit.y)/(firstHit.x-secondHit.x);
  double slope13 = (firstHit.y-thirdHit.y)/(firstHit.x-thirdHit.x);

  if(TMath::Abs(slope12-slope13) < 1e-5) return true;


  return false;
}

bool RingFinder::isPhoton(hit *ahit, int i)
{
  if(ahit->get_pid()->at(i)==0) return true;
  else return false;
}

bool RingFinder::isReflection(hit *ahit, int i)
{
  if(ahit->get_out_pz()->at(i)<0.) return true;
  else return false;
}

bool RingFinder::isOnAerogel(hit *ahit, int i)
{
  // if(ahit->get_out_z()->at(i)>=63.5874 && ahit->get_out_z()->at(i)<=96.5876) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 1) return true;
  else return false;
}

bool RingFinder::isOnPhotonSensor(hit *ahit, int i)
{
  // double out_z = ahit->get_out_z()->at(i);
  // if(out_z > 253.6624+3.0 && out_z < 255.1626+3.0) return true;
  // else return false;
  const int detector_id = ahit->get_id()->at(i);
  if(detector_id == 2) return true;
  else return false;
}

double RingFinder::GausSmearing(TF1 *f_gaus)
{
  double delta_pos = f_gaus->GetRandom();
  return delta_pos;
}

bool RingFinder::isInSensorPlane(double out_x, double out_y)
{
  if( !(TMath::Abs(out_x) > 2.5 && TMath::Abs(out_x) < mRICH::mHalfWidth-2.0) ) return false;
  if( !(TMath::Abs(out_y) > 2.5 && TMath::Abs(out_y) < mRICH::mHalfWidth-2.0) ) return false;
  return true;
}

//--------------------------------------------------------

int RingFinder::Finish()
{
  cout<<endl;
  cout<<"RingFinder::Finish() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(File_mOutPut != NULL){
    File_mOutPut->cd();
    writeHistoMap();
    writeHistoQA();
    File_mOutPut->Close();
  }
  return 0;
}

int RingFinder::writeHistoMap()
{
  h_mPhotonGenerated->Write();

  h_mNumOfPhotons->Write();
  h_mPhotonDist->Write();

  return 0;
}

int RingFinder::writeHistoQA()
{
  h_mXGausSmearing->Write();
  h_mYGausSmearing->Write();
}

//--------------------------------------------------------

int main(int argc, char **argv)
{
  if(argc!=2) return 0;

  const char *input = argv[1];
  string numoflist(input);
  
  string date = "BeamTest_PlusMinus";
  
  cout << "numoflist = " << numoflist.c_str() << endl;
  RingFinder *ringfinder = new RingFinder(numoflist,date);
  
  ringfinder->Init();
  ringfinder->Make();
  ringfinder->Finish();

  cout << "This is the end of RingFinder!!!" << endl;
  
  return 0;
}

