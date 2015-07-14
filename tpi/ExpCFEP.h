#ifndef _TPF_ExpCFEP_
#define _TPF_ExpCFEP_

#include <TH2D.h>
#include <TMath.h>

class ExpCFEP
{
 public:
  ExpCFEP(const char *title, int bineta, int binphi);
  ~ExpCFEP();
  
  void AddRealPair(double deta, double dphi, double weight);
  void AddMixedPair(double deta, double dphi, double weight);
  
  void Write();
  void SetOff(Bool_t dooff);

 private:

  void MakeHistos(const char *title, int bineta, int binphi);

  int isoff;

  TH2D *cnump;
  TH2D *cdenp;
};

#endif
