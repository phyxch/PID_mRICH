#include "../include/Utility.h"
#include "TMath.h"

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

int Utility::get_indexSpaceX(float vx)
{
  // const float vx_start = -57.0; // aerogel halfx is 55.25
  // const float vx_stop = 57.0; // this is to make sure 0 is in cente of [-1,1] and minimize edge effect
  const float delta_vx = 5.0; // in mm
  const float vx_start = -5.0;
  const float vx_stop = 5.0;
  const int NumOfIndex = (int)(vx_stop - vx_start)/delta_vx;

  int index = -1;
  if(vx > vx_start && vx < vx_stop)
  {
    for(int i_index = 0; i_index < NumOfIndex; ++i_index)
    {
      if(vx > vx_start+i_index*delta_vx && vx <= vx_start+(i_index+1)*delta_vx)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexSpaceY(float vy)
{
  // const float vy_start = -57.0; // aerogel halfy is 55.25
  // const float vy_stop = 57.0;
  const float delta_vy = 5.0;
  const float vy_start = -5.0; // aerogel halfy is 55.25
  const float vy_stop = 5.0;
  const int NumOfIndex = (int)(vy_stop - vy_start)/delta_vy;

  int index = -1;
  if(vy > vy_start && vy < vy_stop)
  {
    for(int i_index = 0; i_index < NumOfIndex; ++i_index)
    {
      if(vy > vy_start+i_index*delta_vy && vy <= vy_start+(i_index+1)*delta_vy)
      {
	index = i_index;
      }
    }
  }

  return index;
}

int Utility::get_indexMomentumP(float px, float py, float pz)
{
  const float momentum = TMath::Sqrt(px*px+py*py+pz*pz); // in GeV
  // const float p_start = 2.5;
  // const float p_stop  = 15.5;
  const float p_start = 4.0;
  const float p_stop  = 6.0;
  const float delta_p = 0.2;
  const int NumOfIndex = (int)(p_stop - p_start)/delta_p;

  int index = -1;
  if(momentum > p_start && momentum < p_stop)
  {
    for(int i_index = 0; i_index < NumOfIndex; ++i_index)
    {
      if(momentum > p_start+i_index*delta_p && momentum <= p_start+(i_index+1)*delta_p)
      {
	index = i_index;
      }
    }
  }

  return index;
}
