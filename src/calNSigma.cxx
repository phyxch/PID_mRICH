
#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <math.h> 
#include "string.h"
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include "../include/calNSigma.h"
#include "../include/Utility.h"
#include "../include/mRICH.h"

using namespace std;

calNSigma::calNSigma(string date, string outputfile)
{
  cout<<"calNSigma::calNSigma() ----- Constructor ! ------"<<endl;
  mDate = date;
  mOutPutFile = outputfile;
  utility = new Utility();
}

calNSigma::~calNSigma()
{
  cout<<"calNSigma::~calNSigma() ----- Release memory ! ------"<<endl;
  delete utility;
  delete mFile_OutPut;
}

int calNSigma::Init()
{
  cout<<"calNSigma::init(), create output file: "<< mOutPutFile.c_str() <<endl;
  mFile_OutPut = new TFile(mOutPutFile.c_str(),"RECREATE");

  mInPutFile = Form("/work/eic/xusun/output/probability/PID_prob_%s.root",mDate.c_str());
  initHistoMap_Likelihood();
  initHistoMap_NSigma();

  return 0;
}

int calNSigma::initHistoMap_Likelihood()
{
  cout << "calNSigma::initHistoMap(), initialize histogram for likelihood difference: " << endl;

  cout<<"calLikelihood::init(), read likelihood difference file: "<< mInPutFile.c_str() <<endl;
  mFile_InPut = TFile::Open(mInPutFile.c_str());

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
	    h_mLikelihoodDiff[key_likelihood_first] = (TH2D*)mFile_InPut->Get(key_likelihood_first.c_str());

	    string key_likelihood_second = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_likelihood_second.c_str() << endl;
	    h_mLikelihoodDiff[key_likelihood_second] = (TH2D*)mFile_InPut->Get(key_likelihood_second.c_str());
	  }
	}
      }
    }
  }

  return 0;
}

int calNSigma::initHistoMap_NSigma()
{
  cout << "calNSigma::initHistoMap(), initialize histogram for PID probability: " << endl;
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
	    string key_nsigma_first = utility->gen_KeyNSigma(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    // cout << key_nsigma_first.c_str() << endl;
	    h_mNSigma[key_nsigma_first] = new TH1D(key_nsigma_first.c_str(),key_nsigma_first.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);

	    string key_nsigma_second = utility->gen_KeyNSigma(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_nsigma_second.c_str() << endl;
	    h_mNSigma[key_nsigma_second] = new TH1D(key_nsigma_second.c_str(),key_nsigma_second.c_str(),mRICH::mNumOfIndexMomentumP,mRICH::mMomP_start,mRICH::mMomP_stop);
	  }
	}
      }
    }
  }

  return 0;
}

int calNSigma::writeHistoMap_NSigma()
{
  cout << "calNSigma::writeHistoMap(), write histogram for PID probability: " << endl;
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
	    string key_nsigma_first = utility->gen_KeyNSigma(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,1);
	    // cout << key_nsigma_first.c_str() << endl;
	    h_mNSigma[key_nsigma_first]->Write();

	    string key_nsigma_second = utility->gen_KeyNSigma(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,2);
	    // cout << key_nsigma_second.c_str() << endl;
	    h_mNSigma[key_nsigma_second]->Write();
	  }
	}
      }
    }
  }

  return 0;
}

int calNSigma::Make()
{
  cout << "calNSigma:Make calculate NSigma separation power: " << endl;
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
	    for(int i_pt =0; i_pt < mRICH::mNumOfIndexMomentumP; ++i_pt)
	    {
	      for(int i_rank = 0; i_rank < 2; ++i_rank)
	      {
		string key_likelihood = utility->gen_KeyLikelihood(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,i_rank+1);
		TH1D *h_likelihood_projY = (TH1D*)h_mLikelihoodDiff[key_likelihood]->ProjectionY("h_likelihood_projY",i_pt+1,i_pt+1);
		if(h_likelihood_projY->GetEntries() > 0)
		{
		  double mean  = h_likelihood_projY->GetMean();
		  double width = h_likelihood_projY->GetStdDev();

		  TF1 *f_gaus = new TF1("f_gaus","gaus",-500,500);
		  f_gaus->SetParameter(0,100);
		  f_gaus->SetParameter(1,mean);
		  f_gaus->SetParameter(2,width);
		  f_gaus->SetRange(mean-5.0*width,mean+5.0*width);
		  h_likelihood_projY->Fit(f_gaus,"NR");

		  double mean_diff = f_gaus->GetParameter(1);
		  double width_diff = f_gaus->GetParameter(2);

		  // double sigma_diff = TMath::Sqrt(2.0*mean_diff);
		  // double err_diff = width_diff/TMath::Sqrt(2.0*mean_diff);
		  double sigma_diff = TMath::Sqrt(2.0*mean);
		  double err_diff = width/TMath::Sqrt(2.0*mean);

		  string key_nsigma = utility->gen_KeyNSigma(mRICH::mPIDArray[i_pid],i_vx,i_vy,i_theta,i_phi,i_rank+1);
		  h_mNSigma[key_nsigma]->SetBinContent(i_pt+1,sigma_diff);
		  h_mNSigma[key_nsigma]->SetBinError(i_pt+1,err_diff);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

int calNSigma::Finish()
{
  cout<<endl;
  cout<<"calNSigma::end() ----- Write out tree and histogram to files !------"<<endl;
  cout<<"This is the end of this program !"<<endl;
  if(mFile_OutPut!= NULL){
    mFile_OutPut->cd();
    writeHistoMap_NSigma();
    mFile_OutPut->Close();
  }
  return 0;
}

// This is the main function 
int main()
{
  string date = "May23_2018";
  string outputfile = Form("/work/eic/xusun/output/probability/PID_nSigma_%s.root",date.c_str());

  calNSigma *mcalNSigma = new calNSigma(date,outputfile);
  
  mcalNSigma->Init();
  mcalNSigma->Make();
  mcalNSigma->Finish();

  cout << "This is the end of pid!!!" << endl;
  
  return 0;
}
