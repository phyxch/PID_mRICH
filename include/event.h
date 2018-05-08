//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the event information 
// for EIC RICH detector


#ifndef EVENT_H
#define EVENT_H


#include <vector>
#include <TChain.h>

using namespace std;

class event
{
 public:
  event() {}
  event(TChain *);
  virtual ~event() {}
  
  vector<int>     *get_pid(){ return pid;}
  vector<double>  *get_px(){ return px;}
  vector<double>  *get_py(){ return py;}
  vector<double>  *get_pz(){ return pz;}
  vector<double>  *get_vx(){ return vx;}
  vector<double>  *get_vy(){ return vy;}
  vector<double>  *get_vz(){ return vz;}
  
 private:
  vector<int>     *pid;
  vector<double>  *px;
  vector<double>  *py;
  vector<double>  *pz;
  vector<double>  *vx;
  vector<double>  *vy;
  vector<double>  *vz;
};
#endif
