//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the hit information 
// for EIC RICH detector


#include "../include/hit.h"

hit::hit(TChain *fchain)
{
  id=0;
  hitn=0;
  pid=0;
  mpid=0;
  tid=0;
  mtid=0;
  otid=0;
  trackE=0;
  totEdep=0;
  avg_x=0;
  avg_y=0;
  avg_z=0;
  avg_lx=0;
  avg_ly=0;
  avg_lz=0;
  px=0;
  py=0;
  pz=0;
  vx=0;
  vy=0;
  vz=0;
  mvx=0;
  mvy=0;
  mvz=0;
  avg_t=0;
  in_x=0;
  in_y=0;
  in_z=0;
  in_px=0;
  in_py=0;
  in_pz=0;
  in_t=0;
  out_x=0;
  out_y=0;
  out_z=0;
  out_px=0;
  out_py=0;
  out_pz=0;
  out_t=0;
  
  fchain->SetMakeClass(1);
  fchain->SetBranchAddress("id", &id);
  fchain->SetBranchAddress("hitn", &hitn);
  fchain->SetBranchAddress("pid", &pid);
  fchain->SetBranchAddress("mpid", &mpid);
  fchain->SetBranchAddress("tid", &tid);
  fchain->SetBranchAddress("mtid", &mtid);
  fchain->SetBranchAddress("otid", &otid);
  fchain->SetBranchAddress("trackE", &trackE);
  fchain->SetBranchAddress("totEdep", &totEdep);
  fchain->SetBranchAddress("avg_x", &avg_x);
  fchain->SetBranchAddress("avg_y", &avg_y);
  fchain->SetBranchAddress("avg_z", &avg_z);
  fchain->SetBranchAddress("avg_lx", &avg_lx);
  fchain->SetBranchAddress("avg_ly", &avg_ly);
  fchain->SetBranchAddress("avg_lz", &avg_lz);
  fchain->SetBranchAddress("px", &px);
  fchain->SetBranchAddress("py", &py);
  fchain->SetBranchAddress("pz", &pz);
  fchain->SetBranchAddress("vx", &vx);
  fchain->SetBranchAddress("vy", &vy);
  fchain->SetBranchAddress("vz", &vz);
  fchain->SetBranchAddress("mvx", &mvx);
  fchain->SetBranchAddress("mvy", &mvy);
  fchain->SetBranchAddress("mvz", &mvz);
  fchain->SetBranchAddress("avg_t", &avg_t);
  fchain->SetBranchAddress("in_x", &in_x);
  fchain->SetBranchAddress("in_y", &in_y);
  fchain->SetBranchAddress("in_z", &in_z);
  fchain->SetBranchAddress("in_px", &in_px);
  fchain->SetBranchAddress("in_py", &in_py);
  fchain->SetBranchAddress("in_pz", &in_pz);
  fchain->SetBranchAddress("in_t", &in_t);
  fchain->SetBranchAddress("out_x", &out_x);
  fchain->SetBranchAddress("out_y", &out_y);
  fchain->SetBranchAddress("out_z", &out_z);
  fchain->SetBranchAddress("out_px", &out_px);
  fchain->SetBranchAddress("out_py", &out_py);
  fchain->SetBranchAddress("out_pz", &out_pz);
  fchain->SetBranchAddress("out_t", &out_t);
}
