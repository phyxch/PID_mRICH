#include <string>
#include <map>
#include <vector>
#include <utility>
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

typedef std::map<std::string,TH1D*> TH1DMap;
typedef std::map<std::string,TH2D*> TH2DMap;
typedef std::map<std::string,TH3D*> TH3DMap;
typedef std::map<std::string,TProfile*> TProMap;
typedef std::map<std::string,TGraphAsymmErrors*> TGraMap;
typedef std::map<std::string,TF1*> TF1Map;

typedef std::map<std::string,TH1F*> TH1FMap;
typedef std::map<std::string,TH2F*> TH2FMap;

typedef std::map<std::string,float > floatMap;
typedef std::map<std::string,double> doubleMap;

typedef std::map<int,std::string> strMap;
typedef std::map<int,std::pair<std::string,std::string> > pairMap;

struct hitPosition
{
  double x;
  double y;
}

typedef std::vector<int> intVec;
