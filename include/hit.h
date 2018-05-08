//////*********************************//////
//////***** Author: Dr. Liang Xue *****//////
//////******  Date: 01/06/2015  *******//////
//////*********************************//////

// This class is to define the hit information 
// for EIC RICH detector


#ifndef HIT_H
#define HIT_H


#include <vector>
#include <TChain.h>

using namespace std;

class hit
{
 public:
  hit() {}
  hit(TChain *); 
  virtual ~hit() {}
  
  vector<double>  *get_hitn()   { return hitn;   }
  vector<double>  *get_id()     { return id;     }
  vector<double>  *get_pid()    { return pid;    }
  vector<double>  *get_mpid()   {return mpid;    }
  vector<double>  *get_tid()    {return tid;     }
  vector<double>  *get_mtid()   { return mtid;   }
  vector<double>  *get_otid()   { return otid;   }
  vector<double>  *get_trackE() { return trackE; }
  vector<double>  *get_totEdep(){ return totEdep;}
  vector<double>  *get_avg_x()  { return avg_x;  }
  vector<double>  *get_avg_y()  { return avg_y;  }
  vector<double>  *get_avg_z()  { return avg_z;  }
  vector<double>  *get_avg_lx() { return avg_lx; }
  vector<double>  *get_avg_ly() { return avg_ly; }
  vector<double>  *get_avg_lz() { return avg_lz; }
  vector<double>  *get_px()     { return px;     }
  vector<double>  *get_py()     { return py;     }
  vector<double>  *get_pz()     { return pz;     }
  vector<double>  *get_vx()     { return vx;     }
  vector<double>  *get_vy()     { return vy;     }
  vector<double>  *get_vz()     { return vz;     }
  vector<double>  *get_mvx()    { return mvx;    }
  vector<double>  *get_mvy()    { return mvy;    }
  vector<double>  *get_mvz()    { return mvz;    }
  vector<double>  *get_avg_t()  { return avg_t;  }
  vector<double>  *get_in_x()   { return in_x;   }
  vector<double>  *get_in_y()   { return in_y;   }
  vector<double>  *get_in_z()   { return in_z;   }
  vector<double>  *get_in_px()  { return in_px;  }
  vector<double>  *get_in_py()  { return in_py;  }
  vector<double>  *get_in_pz()  { return in_pz;  }
  vector<double>  *get_in_t()   { return in_t;   }
  vector<double>  *get_out_x()  { return out_x;  }
  vector<double>  *get_out_y()  { return out_y;  }
  vector<double>  *get_out_z()  { return out_z;  }
  vector<double>  *get_out_px() { return out_px; }
  vector<double>  *get_out_py() { return out_py; }
  vector<double>  *get_out_pz() { return out_pz; }
  vector<double>  *get_out_t()  { return out_t;  }
  
 private:
  vector<double>  *id;
  vector<double>  *hitn;
  vector<double>  *pid;
  vector<double>  *mpid;
  vector<double>  *tid;
  vector<double>  *mtid;
  vector<double>  *otid;
  vector<double>  *trackE;
  vector<double>  *totEdep;
  vector<double>  *avg_x;
  vector<double>  *avg_y;
  vector<double>  *avg_z;
  vector<double>  *avg_lx;
  vector<double>  *avg_ly;
  vector<double>  *avg_lz;
  vector<double>  *px;
  vector<double>  *py;
  vector<double>  *pz;
  vector<double>  *vx;
  vector<double>  *vy;
  vector<double>  *vz;
  vector<double>  *mvx;
  vector<double>  *mvy;
  vector<double>  *mvz;
  vector<double>  *avg_t;
  vector<double>  *in_x;
  vector<double>  *in_y;
  vector<double>  *in_z;
  vector<double>  *in_px;
  vector<double>  *in_py;
  vector<double>  *in_pz;
  vector<double>  *in_t;
  vector<double>  *out_x;
  vector<double>  *out_y;
  vector<double>  *out_z;
  vector<double>  *out_px;
  vector<double>  *out_py;
  vector<double>  *out_pz;
  vector<double>  *out_t;
};
#endif
