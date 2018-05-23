#include "../include/Utility.h"
#include "../include/mRICH.h"
#include "TMath.h"
#include "TString.h"

Utility::Utility()
{
  creatParticleIdMap();
  creatParticleMisIdMap();
}

Utility::~Utility()
{
  /* */
}

void Utility::creatParticleIdMap()
{
  mParticleIdMap.clear();
  mParticleIdMap[211]   = "piplus";
  mParticleIdMap[321]   = "Kplus";
  mParticleIdMap[2212]  = "proton";
  mParticleIdMap[-211]  = "piminus";
  mParticleIdMap[-321]  = "Kminus";
  mParticleIdMap[-2212] = "antiproton";
}

void Utility::creatParticleMisIdMap()
{
  mParticleMisIdMap.clear();
  mParticleMisIdMap[211]    = std::pair<std::string,std::string>("Kplus","proton");
  mParticleMisIdMap[321]    = std::pair<std::string,std::string>("piplus","proton");
  mParticleMisIdMap[2212]   = std::pair<std::string,std::string>("piplus","Kplus");
  mParticleMisIdMap[-211]   = std::pair<std::string,std::string>("Kminus","antiproton");
  mParticleMisIdMap[-321]   = std::pair<std::string,std::string>("piminus","antiproton");
  mParticleMisIdMap[-2212]  = std::pair<std::string,std::string>("piminus","Kminus");
}

std::string Utility::get_IdentifiedParticle(int pid)
{
  return mParticleIdMap[pid];
}

std::pair<std::string,std::string> Utility::get_misIdentifiedParticle(int pid)
{
  return mParticleMisIdMap[pid];
}

int Utility::get_indexSpaceX(double vx)
{
  // const double delta_vx = (mRICH::mVx_stop-mRICH::mVx_start)/mRICH::mNumOfIndexSpaceX; 

  int index = -1;
  if(vx >= mRICH::mVx_start && vx < mRICH::mVx_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexSpaceX; ++i_index)
    {
      // if(vx >= mRICH::mVx_start+i_index*delta_vx && vx < mRICH::mVx_start+(i_index+1)*delta_vx)
      if(vx >= mRICH::mBin_Vx[i_index]-mRICH::mDelta_Vx && vx < mRICH::mBin_Vx[i_index]+mRICH::mDelta_Vx)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexSpaceY(double vy)
{
  // const double delta_vy = (mRICH::mVy_stop-mRICH::mVy_start)/mRICH::mNumOfIndexSpaceY; 

  int index = -1;
  if(vy >= mRICH::mVy_start && vy < mRICH::mVy_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexSpaceY; ++i_index)
    {
      // if(vy >= mRICH::mVy_start+i_index*delta_vy && vy < mRICH::mVy_start+(i_index+1)*delta_vy)
      if(vy >= mRICH::mBin_Vy[i_index]-mRICH::mDelta_Vy && vy < mRICH::mBin_Vy[i_index]+mRICH::mDelta_Vy)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexMomentumP(double px, double py, double pz)
{
  const double momentum = TMath::Sqrt(px*px+py*py+pz*pz); // in GeV

  // const double delta_p = (mRICH::mMomP_stop-mRICH::mMomP_start)/mRICH::mNumOfIndexMomentumP; 

  int index = -1;
  if(momentum >= mRICH::mMomP_start && momentum < mRICH::mMomP_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexMomentumP; ++i_index)
    {
      // if(momentum >= mRICH::mMomP_start+i_index*delta_p && momentum < mRICH::mMomP_start+(i_index+1)*delta_p)
      if(momentum >= mRICH::mBin_MomP[i_index]-mRICH::mDelta_MomP && momentum < mRICH::mBin_MomP[i_index]+mRICH::mDelta_MomP)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexMomentumTheta(double px, double py, double pz)
{
  const double momentum = TMath::Sqrt(px*px+py*py+pz*pz); // in GeV

  // const double delta_theta = (mRICH::mMomTheta_stop-mRICH::mMomTheta_start)/mRICH::mNumOfIndexMomentumTheta; 
  const double theta = TMath::ACos(pz/momentum)*mRICH::DEG;    //in deg

  int index = -1;
  if(theta >= mRICH::mMomTheta_start && theta < mRICH::mMomTheta_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexMomentumTheta; ++i_index)
    {
      // if(theta >= mRICH::mMomTheta_start+i_index*delta_theta && theta < mRICH::mMomTheta_start+(i_index+1)*delta_theta)
      if(theta >= mRICH::mBin_MomTheta[i_index]-mRICH::mDelta_MomTheta && theta < mRICH::mBin_MomTheta[i_index]+mRICH::mDelta_MomTheta)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexMomentumPhi(double px, double py)
{
  // const double delta_phi = (mRICH::mMomPhi_stop-mRICH::mMomPhi_start)/mRICH::mNumOfIndexMomentumPhi; 
  const double phi = TMath::ATan2(py,px)*mRICH::DEG;    //in deg            

  int index = -1;
  if(phi >= mRICH::mMomPhi_start && phi < mRICH::mMomPhi_stop)
  {
    for(int i_index = 0; i_index < mRICH::mNumOfIndexMomentumPhi; ++i_index)
    {
      // if(phi >= mRICH::mMomPhi_start+i_index*delta_phi && phi < mRICH::mMomPhi_start+(i_index+1)*delta_phi)
      if(phi >= mRICH::mBin_MomPhi[i_index]-mRICH::mDelta_MomPhi && phi < mRICH::mBin_MomPhi[i_index]+mRICH::mDelta_MomPhi)
      {
	index = i_index;
      }
    }
  }

  return index;
}


std::string Utility::gen_KeyNumOfEvents(int pid, int index_vx, int index_vy, int index_p, int index_theta, int index_phi)
{
  std::string identifiedParticle = this->get_IdentifiedParticle(pid);
  std::string key_events = Form("h_mNumOfEvents_%s_vx_%d_vy_%d_mom_%d_theta_%d_phi_%d",identifiedParticle.c_str(),index_vx,index_vy,index_p,index_theta,index_phi);

  return key_events;
}

std::string Utility::gen_KeyMassHypo(int pid, int index_vx, int index_vy, int index_p, int index_theta, int index_phi)
{
  std::string identifiedParticle = this->get_IdentifiedParticle(pid);
  std::string key_photon = Form("h_mPhotonDist_%s_vx_%d_vy_%d_mom_%d_theta_%d_phi_%d",identifiedParticle.c_str(),index_vx,index_vy,index_p,index_theta,index_phi);

  return key_photon;
}

std::string Utility::gen_KeyLikelihood(int pid, int index_vx, int index_vy, int index_theta, int index_phi, int rank)
{
  std::string identifiedParticle = this->get_IdentifiedParticle(pid);
  std::pair<std::string,std::string> misIdentifiedParticle = this->get_misIdentifiedParticle(pid);
  std::string key_likelihood = "undifined";

  if(rank == 1)
  {
    key_likelihood = Form("h_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  }
  if(rank == 2)
  {
    key_likelihood = Form("h_mLikelihoodDiff_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  }

  return key_likelihood;
}

std::string Utility::gen_KeyProb(int pid, int index_vx, int index_vy, int index_theta, int index_phi, int rank)
{
  std::string identifiedParticle = this->get_IdentifiedParticle(pid);
  std::pair<std::string,std::string> misIdentifiedParticle = this->get_misIdentifiedParticle(pid);
  std::string key_prob = "undifined";

  if(rank == 0)
  {
    key_prob = Form("h_mProbability_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),identifiedParticle.c_str(),index_vx,index_vy,index_theta,index_phi);
  }
  if(rank == 1)
  {
    key_prob = Form("h_mProbability_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),misIdentifiedParticle.first.c_str(),index_vx,index_vy,index_theta,index_phi);
  }
  if(rank == 2)
  {
    key_prob = Form("h_mProbability_%s_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),misIdentifiedParticle.second.c_str(),index_vx,index_vy,index_theta,index_phi);
  }

  return key_prob;
}

std::string Utility::gen_KeySumOfPID(int pid, int index_vx, int index_vy, int index_theta, int index_phi)
{
  std::string identifiedParticle = this->get_IdentifiedParticle(pid);

  std::string key_sumofpid = Form("h_mSumOfPID_%s_vx_%d_vy_%d_theta_%d_phi_%d",identifiedParticle.c_str(),index_vx,index_vy,index_theta,index_phi);

  return key_sumofpid;
}
