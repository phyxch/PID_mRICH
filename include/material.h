//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the photonsensor
// materials for EIC RICH detector


#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

using namespace std;

class material
{
 public:
  material();
  virtual ~material() {}
  
  double extrapQE_SbKCs(double wl);
  double extrapQE_GaAsP(double wl);
  double extrapQE_GaAs(double wl);
  double extrapQE(double wl);
  
 private:
  vector< pair<double, double> > SbKCs;
  vector< pair<double, double> > GaAsP;
  vector< pair<double, double> > GaAs;
  vector< pair<double, double> > QE_Ping;
};
#endif
