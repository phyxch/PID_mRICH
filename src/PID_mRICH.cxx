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
#include "../include/PID_mRICH.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;

PID_mRICH::PID_mRICH(string date, string outputfile)
{
 cout<<"PID_mRICH::PID_mRICH() ----- Constructor ! ------"<<endl;
 mDate = date;
 mOutPutFile = outputfile;
 utility = new Utility();
}

PID_mRICH::~PID_mRICH()
{
 cout<<"PID_mRICH::~PID_mRICH() ----- Release memory ! ------"<<endl;
 delete utility;
 delete mFile_OutPut;
 delete mChainInPut;
}

int PID_mRICH::Init()
{
 cout<<"PID_mRICH::init(), create output file: "<< mOutPutFile.c_str() <<endl;
 mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");
 initChain();
 initHistoMap_Likelihood();
 initHistoMap_Probability();

  return 0;
}

int PID_mRICH::initChain()
{
 string inputdir = "/work/eic/xusun/output/likelihood/";
 string InPutList = Form("/work/eic/xusun/list/probability/mRICH_Prob_%s.list",mDate.c_str());

 mChainInPut = new TChain("LLTreeDst");

 if (!InPutList.empty())   // if input file is ok
 {
   cout << "Open input probability file list" << endl;
   ifstream in(InPutList.c_str());  // input stream
   if(in)
   {
     cout << "input file probability list is ok" << endl;
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
	 mChainInPut->AddFile(addfile.c_str(),-1,"LLTreeDst");
	 Long64_t file_entries = mChainInPut->GetEntries();
	 cout << "File added to data chain: " << addfile.c_str() << " with " << (file_entries-entries_save) << " entries" << endl;
	 entries_save = file_entries;
       }
     }
   }
   else
   {
     cout << "WARNING: input probability file input is problemtic" << endl;
   }
 }

 mChainInPut->SetBranchAddress("pid", &mPid);
 mChainInPut->SetBranchAddress("px", &mPx);
 mChainInPut->SetBranchAddress("py", &mPy);
 mChainInPut->SetBranchAddress("pz", &mPz);
 mChainInPut->SetBranchAddress("vx", &mVx);
 mChainInPut->SetBranchAddress("vy", &mVy);
 mChainInPut->SetBranchAddress("vz", &mVz);
 mChainInPut->SetBranchAddress("theta", &mTheta);
 mChainInPut->SetBranchAddress("phi", &mPhi);
 mChainInPut->SetBranchAddress("nhit", &mNHit);
 mChainInPut->SetBranchAddress("nhitAegl", &mNHitAegl);
 mChainInPut->SetBranchAddress("nhitPhoDet", &mNHitPhoDet);
 mChainInPut->SetBranchAddress("nelSbKCs", &mNelSbKCs);
 mChainInPut->SetBranchAddress("nelGaAsP", &mNelGaAsP);
 mChainInPut->SetBranchAddress("nelGaAs", &mNelGaAs);
 mChainInPut->SetBranchAddress("Lpion", &mLpion);
 mChainInPut->SetBranchAddress("LKaon", &mLKaon);
 mChainInPut->SetBranchAddress("Lproton", &mLproton);

 long NumOfEvents = (long)mChainInPut->GetEntries();
 cout << "total number of events: " << NumOfEvents << endl;

 return 0;
}

int PID_mRICH::initHistoMap_Likelihood()
{
  cout << "PID_mRICH::initHistoMap(), initialize histogram for likelihood difference: " << endl;
  // initialization of h_mLikelihoodDiff
  for(int i_pid = 0; i_pid < mRICH::mNumOfParType; ++i_pid)
  {
    for(int i_vx = 0; i_vx < mRICH::mNumOfIndexSpaceX; ++i_vx)
    {
      for(int i_vy = 0; i_vy < mRICH::mNumOfIndexSpaceY; ++i_vy)
      {
	for(int i_theta = 0; i_theta < mRICH::mNumOfIndexMomentumTheta; ++i_theta)
	{
	  for(int i_phi = 0; i_phi < mRICH::mNumOfIndexMomentumPhi; ++i_phi)
	  {
	    string key_likelihood_first = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    // cout << key_likelihood_first.c_str() << endl;
	    h_mLikelihoodDiff[key_likelihood_first] = new TH2D(key_likelihood_first.c_str(),key_likelihood_first.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop,mRICH::mNumOfLieklihood,mRICH::mLieklihood_start,mRICH::mLieklihood_stop);

	    string key_likelihood_second = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_likelihood_second.c_str() << endl;
	    h_mLikelihoodDiff[key_likelihood_second] = new TH2D(key_likelihood_second.c_str(),key_likelihood_second.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop,mRICH::mNumOfLieklihood,mRICH::mLieklihood_start,mRICH::mLieklihood_stop);
	  }
	}
      }
    }
  }

  return 0;
}

int PID_mRICH::writeHistoMap_Likelihood()
{
  cout << "PID_mRICH::writeHistoMap(), write histogram for likelihood difference: " << endl;
  // initialization of h_mLikelihoodDiff
  for(int i_pid = 0; i_pid < mRICH::mNumOfParType; ++i_pid)
  {
    for(int i_vx = 0; i_vx < mRICH::mNumOfIndexSpaceX; ++i_vx)
    {
      for(int i_vy = 0; i_vy < mRICH::mNumOfIndexSpaceY; ++i_vy)
      {
	for(int i_theta = 0; i_theta < mRICH::mNumOfIndexMomentumTheta; ++i_theta)
	{
	  for(int i_phi = 0; i_phi < mRICH::mNumOfIndexMomentumPhi; ++i_phi)
	  {
	    string key_likelihood_first = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    h_mLikelihoodDiff[key_likelihood_first]->Write();

	    string key_likelihood_second = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    h_mLikelihoodDiff[key_likelihood_second]->Write();
	  }
	}
      }
    }
  }

  return 0;
}

int PID_mRICH::initHistoMap_Probability()
{
  cout << "PID_mRICH::initHistoMap(), initialize histogram for PID probability: " << endl;
  // initialization of h_mProbability
  for(int i_pid = 0; i_pid < mRICH::mNumOfParType; ++i_pid)
  {
    for(int i_vx = 0; i_vx < mRICH::mNumOfIndexSpaceX; ++i_vx)
    {
      for(int i_vy = 0; i_vy < mRICH::mNumOfIndexSpaceY; ++i_vy)
      {
	for(int i_theta = 0; i_theta < mRICH::mNumOfIndexMomentumTheta; ++i_theta)
	{
	  for(int i_phi = 0; i_phi < mRICH::mNumOfIndexMomentumPhi; ++i_phi)
	  {
	    string key_prob_identified = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,0);
	    // cout << key_prob_identified.c_str() << endl;
	    h_mProbability[key_prob_identified] = new TH1D(key_prob_identified.c_str(),key_prob_identified.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
	    h_mProbability[key_prob_identified]->Sumw2();

	    string key_prob_misIdentified_first = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    // cout << key_prob_misIdentified_first.c_str() << endl;
	    h_mProbability[key_prob_misIdentified_first] = new TH1D(key_prob_misIdentified_first.c_str(),key_prob_misIdentified_first.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
	    h_mProbability[key_prob_misIdentified_first]->Sumw2();

	    string key_prob_misIdentified_second = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_prob_misIdentified_second.c_str() << endl;
	    h_mProbability[key_prob_misIdentified_second] = new TH1D(key_prob_misIdentified_second.c_str(),key_prob_misIdentified_second.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
	    h_mProbability[key_prob_misIdentified_second]->Sumw2();

	    // total number of particles to be identified
	    string key_sumofpid = utility->gen_KeySumOfPID(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi);
	    h_mSumOfPID[key_sumofpid] = new TH1D(key_sumofpid.c_str(),key_sumofpid.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
	    h_mSumOfPID[key_sumofpid]->Sumw2();
	  }
	}
      }
    }
  }

  return 0;
}

int PID_mRICH::writeHistoMap_Probability()
{
  cout << "PID_mRICH::writeHistoMap(), write histogram for PID probability: " << endl;
  // initialization of h_mProbability
  for(int i_pid = 0; i_pid < mRICH::mNumOfParType; ++i_pid)
  {
    for(int i_vx = 0; i_vx < mRICH::mNumOfIndexSpaceX; ++i_vx)
    {
      for(int i_vy = 0; i_vy < mRICH::mNumOfIndexSpaceY; ++i_vy)
      {
	for(int i_theta = 0; i_theta < mRICH::mNumOfIndexMomentumTheta; ++i_theta)
	{
	  for(int i_phi = 0; i_phi < mRICH::mNumOfIndexMomentumPhi; ++i_phi)
	  {
	    string key_prob_identified = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,0);
	    // cout << key_prob_identified.c_str() << endl;
	    h_mProbability[key_prob_identified]->Write();

	    string key_prob_misIdentified_first = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    // cout << key_prob_misIdentified_first.c_str() << endl;
	    h_mProbability[key_prob_misIdentified_first]->Write();

	    string key_prob_misIdentified_second = utility->gen_KeyProb(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_prob_misIdentified_second.c_str() << endl;
	    h_mProbability[key_prob_misIdentified_second]->Write();

	    // total number of particles to be identified
	    string key_sumofpid = utility->gen_KeySumOfPID(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi);
	    h_mSumOfPID[key_sumofpid]->Write();
	  }
	}
      }
    }
  }

  return 0;
}

std::pair<double,double> PID_mRICH::get_LikelihoodDiff(int pid, double Lpion, double LKaon, double Lproton)
{
  if(TMath::Abs(pid) ==  211)   return std::pair<double,double>(Lpion-LKaon,Lpion-Lproton);
  if(TMath::Abs(pid) ==  321)   return std::pair<double,double>(LKaon-Lpion,LKaon-Lproton);
  if(TMath::Abs(pid) ==  2212)  return std::pair<double,double>(Lproton-Lpion,Lproton-LKaon);
  else return std::pair<double,double>(-999.9,-999.9);
}


std::pair<int,int> PID_mRICH::get_misPID(int pid)
{
   if(pid == 211)   return std::pair<int,int>(321,2212);
   if(pid == 321)   return std::pair<int,int>(211,2212);
   if(pid == 2212)  return std::pair<int,int>(211,321);
   if(pid == -211)  return std::pair<int,int>(-321,-2212);
   if(pid == -321)  return std::pair<int,int>(-211,-2212);
   if(pid == -2212) return std::pair<int,int>(-211,-321);
}

int PID_mRICH::get_rank(int pid, double Lpion, double LKaon, double Lproton)
{
  int rank = -1;
  std::pair<double,double> likelihood_diff = this->get_LikelihoodDiff(pid,Lpion,LKaon,Lproton);

  std::pair<int,int> misPID = this->get_misPID(pid);
  std::pair<double,double> likelihood_diff_first = this->get_LikelihoodDiff(misPID.first,Lpion,LKaon,Lproton);
  std::pair<double,double> likelihood_diff_second = this->get_LikelihoodDiff(misPID.second,Lpion,LKaon,Lproton);

  if(likelihood_diff.first > 0.0 && likelihood_diff.second > 0.0) rank = 0;
  else if(likelihood_diff_first.first > 0.0 && likelihood_diff_first.second > 0.0) rank = 1;
  else if(likelihood_diff_second.first > 0.0 && likelihood_diff_second.second > 0.0) rank = 2;

  return rank;
}

int PID_mRICH::Make()
{
  long NumOfEvents = (long)mChainInPut->GetEntries();

  mChainInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  for(long i_event = 0; i_event < NumOfEvents; ++i_event)
  // for(long i_event = 0; i_event < 100; ++i_event)
  {
    if (!mChainInPut->GetEntry(i_event)) // take the event -> information is stored in event
      break;

    if(i_event%100==0) cout << "processing event:  " << i_event << " ;"<<endl;

    const int indexSpaceX = utility->get_indexSpaceX(mVx);
    const int indexSpaceY = utility->get_indexSpaceY(mVy);
    const int indexMomentumTheta = utility->get_indexMomentumTheta(mPx,mPy,mPz);
    const int indexMomentumPhi = utility->get_indexMomentumPhi(mPx,mPy);
    const double momentum = TMath::Sqrt(mPx*mPx+mPy*mPy+mPz*mPz); // in GeV

    std::pair<double,double> likelihood_diff = this->get_LikelihoodDiff(mPid,mLpion,mLKaon,mLproton);

    //---- likelihood difference QA ----
    string key_likelihood_first = utility->gen_KeyLikelihood(mPid,indexSpaceX,indexSpaceY,indexMomentumTheta,indexMomentumPhi,1);
    h_mLikelihoodDiff[key_likelihood_first]->Fill(momentum,likelihood_diff.first);

    string key_likelihood_second = utility->gen_KeyLikelihood(mPid,indexSpaceX,indexSpaceY,indexMomentumTheta,indexMomentumPhi,2);
    h_mLikelihoodDiff[key_likelihood_second]->Fill(momentum,likelihood_diff.second);
    //---- likelihood difference QA ----

    int rank = this->get_rank(mPid,mLpion,mLKaon,mLproton);
    if(rank < 0) continue;

    string key_prob = utility->gen_KeyProb(mPid,indexSpaceX,indexSpaceY,indexMomentumTheta,indexMomentumPhi,rank);
    h_mProbability[key_prob]->Fill(momentum);

    string key_sumofpid = utility->gen_KeySumOfPID(mPid,indexSpaceX,indexSpaceY,indexMomentumTheta,indexMomentumPhi);
    h_mSumOfPID[key_sumofpid]->Fill(momentum);
  }

  return 0;
}

int PID_mRICH::Finish()
{
  cout<<endl;
  cout<<"PID_mRICH::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(mFile_OutPut!= NULL){
    mFile_OutPut->cd();
    writeHistoMap_Likelihood();
    writeHistoMap_Probability();
    mFile_OutPut->Close();
  }
  return 0;
}

////// This is the main function 
int main()
{
  string date = "May15_2018";
  string outputfile = Form("/work/eic/xusun/output/probability/PID_prob_%s.root",date.c_str());

  PID_mRICH *mPID_mRICH = new PID_mRICH(date,outputfile);
  
  mPID_mRICH->Init();
  mPID_mRICH->Make();
  mPID_mRICH->Finish();
  
  return 0;
}
