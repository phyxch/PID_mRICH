#ifndef mRICH_h
#define mRICH_h

#include <string>
#include "TMath.h"

namespace mRICH
{
  // vterx segamentation
  // const double mVx_start = -57.0; // aerogel halfx is 55.25
  // const double mVx_stop = 57.0; // this is to make sure 0 is in cente of [-1,1] and minimize edge effect
  const double mVx_start = -5.0;
  const double mVx_stop = 5.0;
  const int mNumOfIndexSpaceX = 2;

  // const double mVy_start = -57.0;
  // const double mVy_stop = 57.0;
  const double mVy_start = -5.0;
  const double mVy_stop = 5.0;
  const int mNumOfIndexSpaceY = 2;

  // momentum segamentation
  // const double p_start = 2.5;
  // const double p_stop  = 15.5;
  const double mMomP_start = 4.0;
  const double mMomP_stop  = 6.0;
  const int mNumOfIndexMomentumP = 10;

  const double DEG=180./TMath::Pi();
  const double mMomTheta_start = 0.0*DEG;
  const double mMomTheta_stop  = 0.25*TMath::Pi()*DEG;
  const int mNumOfIndexMomentumTheta = 1;

  const double mMomPhi_start = -1.0*TMath::Pi()*DEG;
  const double mMomPhi_stop  = TMath::Pi()*DEG;
  const int mNumOfIndexMomentumPhi = 1;

  // initialization for mass hypotheses histogram
  const int mNumOfParType = 6;
  const std::string mPID[6] = {"piplus","Kplus","proton","piminus","Kminus","antiproton"};
  const int mPid[mNumOfParType] = {211,321,2212,-211,-321,-2212};
  const int mNPads=105;   ///// number of photonsenor segmentation pads
  const double mHalfWidth=52.5; // (mm) half width (x) and height (y) of the photon sensor

  const int mLieklihood_start = -50;
  const int mLieklihood_stop  = 450;
  const int mNumOfLieklihood = 500;
}

#endif
