//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the event information 
// for EIC RICH detector


#include "../include/event.h"

event::event(TChain *fchain)
{
  pid = 0;
  px = 0;
  py = 0;
  pz = 0;
  vx = 0;
  vy = 0;
  vz = 0;
  
  fchain->SetMakeClass(1);
  
  fchain->SetBranchAddress("pid", &pid);
  fchain->SetBranchAddress("px", &px);
  fchain->SetBranchAddress("py", &py);
  fchain->SetBranchAddress("pz", &pz);
  fchain->SetBranchAddress("vx", &vx);
  fchain->SetBranchAddress("vy", &vy);
  fchain->SetBranchAddress("vz", &vz);
}
