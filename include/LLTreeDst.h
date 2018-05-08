//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the tree format 


#ifndef LLTREEDST_H
#define LLTREEDST_H

class LLTreeDst
{
 public:
  LLTreeDst() {}
  virtual ~LLTreeDst() {}
  
  int pid;
  float px;
  float py;
  float pz;
  float vx;
  float vy;
  float vz;
  float theta;
  float phi;
  int nhit;
  int nhitAegl;
  int nhitPhoDet;
  int nelSbKCs;
  int nelGaAsP;
  int nelGaAs;
  double Lpion;
  double LKaon;
  double Lproton;
};
#endif
