#ifndef mRICH_h
#define mRICH_h

#include <string>
#include "TMath.h"

namespace mRICH
{
  // vterx segamentation
  const double mVx_start = -57.0; // aerogel halfx is 55.25
  const double mVx_stop = 57.0; // this is to make sure 0 is in cente of [-1,1] and minimize edge effect
  // const double mVx_start = -5.0;
  // const double mVx_stop = 5.0;
  const int mNumOfIndexSpaceX = 1;
  const double mBin_Vx[2] = {0,10};
  // const double mBin_Vx[mNumOfIndexSpaceX] = {-15,-10,0,10,15};
  const double mDelta_Vx = 5;

  const double mVy_start = -57.0;
  const double mVy_stop = 57.0;
  // const double mVy_start = -5.0;
  // const double mVy_stop = 5.0;
  const int mNumOfIndexSpaceY = 1;
  const double mBin_Vy[2] = {0,10};
  // const double mBin_Vy[mNumOfIndexSpaceY] = {-15,-10,0,10,15};
  const double mDelta_Vy = 5;

  // momentum segamentation
  const double mMomP_start = 2.5;
  const double mMomP_stop  = 15.5;
  const int mNumOfIndexMomentumP = 13;
  const double mBin_MomP[mNumOfIndexMomentumP] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
  // const double mMomP_start = 2.5;
  // const double mMomP_stop  = 10.5;
  // const int mNumOfIndexMomentumP = 8;
  // const double mBin_MomP[mNumOfIndexMomentumP] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
  const double mDelta_MomP = 0.5;

  const double DEG=180./TMath::Pi();

  const double mMomTheta_start = 0.0*DEG;
  const double mMomTheta_stop  = 0.25*TMath::Pi()*DEG;
  // const double mMomTheta_start = 0.0*DEG;
  // const double mMomTheta_stop  = TMath::Pi()*DEG/6.0;
  const int mNumOfIndexMomentumTheta = 5;
  const double mBin_MomTheta[mNumOfIndexMomentumTheta] = {2.5,12.5,22.5,32.5,42.5};
  const double mDelta_MomTheta = 2.5;

  const double mMomPhi_start = -1.0*TMath::Pi()*DEG;
  const double mMomPhi_stop  = TMath::Pi()*DEG;
  // const double mMomPhi_start = 0.0*DEG;
  // const double mMomPhi_stop  = TMath::Pi()*DEG/3.0;
  const int mNumOfIndexMomentumPhi = 7;
  const double mBin_MomPhi[mNumOfIndexMomentumPhi] = {-90.0,-60.0,-30.0,0.0,30.0,60.0,90.0};
  const double mDelta_MomPhi = 2.5;

  // initialization for mass hypotheses histogram
  const int mNumOfParType = 6;
  const int mPIDArray[mNumOfParType] = {211,321,2212,-211,-321,-2212};

  // const int mNumOfPixels = 105;   ///// number of photonsenor segmentation pads
  const double mHalfWidth = 52.5; // (mm) half width (x) and height (y) of the glass window
  const int mNumOfPixels = 37; // 16*2 3mm-pixels + 2*2 2mm-glasswindow + 1 1mm-gap
  const double mPixels[mNumOfPixels+1] = {-52.5,-50.5,-47.5,-44.5,-41.5,-38.5,-35.5,-32.5,-29.5,-26.5,-23.5,-20.5,-17.5,-14.5,-11.5,-8.5,-5.5,-2.5,-0.5,0.5,2.5,5.5,8.5,11.5,14.5,17.5,20.5,23.5,26.5,29.5,32.5,35.5,38.5,41.5,44.5,47.5,50.5,52.5};

  const int mLieklihood_start = -50;
  const int mLieklihood_stop  = 450;
  const int mNumOfLieklihood = 500;
}

#endif
