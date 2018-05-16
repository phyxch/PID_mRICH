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
    int get_indexMomentumTheta(double px, double py, double pz);
    int get_indexMomentumPhi(double px, double py);

    // std::string gen_KeyMassHypo(int pid, int index_vx, int index_vy, int index_p, int index_theta, int index_phi);

    // order: 0 for identified | 1 for 1st misIdentified | 2 for 2nd misIdentified
    std::string gen_KeyLikelihood(int pid, int index_vx, int index_vy, int index_theta, int index_phi, int order);
    std::string gen_KeyProb(int pid, int index_vx, int index_vy, int index_theta, int index_phi, int order);


  private:
    strMap mParticleIdMap;
    pairMap mParticleMisIdMap;

    void creatParticleIdMap();
    void creatParticleMisIdMap();
};

#endif
