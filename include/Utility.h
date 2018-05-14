#ifndef Utility_h
#define Utility_h

#include <string>
#include <utility>
#include <map>
#include "type.h"

class Utility
{
  public:
    Utility();
    ~Utility();

    std::string get_IdentifiedParticle(int pid);
    std::pair<std::string,std::string> get_misIdentifiedParticle(int pid);
    int get_indexSpaceX(double vx);
    int get_indexSpaceY(double vy);
    int get_indexMomentumP(double px, double py, double pz);
    // int get_indexMomentumTheta(double px, double py, double pz);
    // int get_indexMomentumPhi(double px, double py, double pz);

  private:
    strMap mParticleIdMap;
    pairMap mParticleMisIdMap;

    void creatParticleIdMap();
    void creatParticleMisIdMap();
};

#endif
