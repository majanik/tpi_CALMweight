#ifndef _TPF_ExpCF1D_
#define _TPF_ExpCF1D_

#include <TH1D.h>

class ExpCF1D
{
 public:
  ExpCF1D(const char *title, int bin, double min, double max);
  ~ExpCF1D();
  
  void AddRealPair(double qv, double weight, int sw=0);
  void AddMixedPair(double qv, double weight, int sw=0);
  
  void Write();
  void SetOff(Bool_t dooff);

 private:

  void MakeHistos(const char *title, int bin, double min, double max);

  int isoff;

  TH1D *cnump;
  TH1D *cdenp;

  TH1D *cnumn;
  TH1D *cdenn;
};

#endif
