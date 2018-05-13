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
    int get_indexSpaceX(float vx);
    int get_indexSpaceY(float vy);
    int get_indexMomentumP(float px, float py, float pz);
    // int get_indexMomentumTheta(float px, float py, float pz);
    // int get_indexMomentumPhi(float px, float py, float pz);

  private:
    strMap mParticleIdMap;
    pairMap mParticleMisIdMap;

    void creatParticleIdMap();
    void creatParticleMisIdMap();
};

#endif
