#ifndef _TPF_ExpCFSH_
#define _TPF_ExpCFSH_

#include "CorrFctnDirectYlm.h"

class ExpCFSH
{
 public:
  ExpCFSH(const char *title, int bin, double min, double max);
  ~ExpCFSH();
  
  void AddRealPair(double qx, double qy, double qz, double weight);
  void AddMixedPair(double qx, double qy, double qz, double weight);
  
  void Write();
  void SetOff(Bool_t dooff);

 private:

  void MakeHistos(const char *title, int bin, double min, double max);

  int isoff;

  CorrFctnDirectYlm *cfun;
};

#endif
