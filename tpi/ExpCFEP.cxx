#include "ExpCFEP.h"

ExpCFEP::ExpCFEP(const char *title, int bineta, int binphi)
{
  isoff = false;
  MakeHistos(title, bineta, binphi);
}

ExpCFEP::~ExpCFEP()
{
  if (cnump) delete cnump;
  if (cdenp) delete cdenp;
}
  
void ExpCFEP::SetOff(Bool_t dooff)
{
  isoff = dooff;
}

void ExpCFEP::MakeHistos(const char *title, int bineta, int binphi)
{
  char htitle[2000];

      double phiL = (-(int)(binphi/4)+0.5)*2.*TMath::Pi()/binphi;
      double phiU = 2*TMath::Pi()+(-(int)(binphi/4)+0.5)*2.*TMath::Pi()/binphi;

  sprintf(htitle, "cnumep%s", title);
  cnump = new TH2D(htitle, htitle, binphi, phiL, phiU, bineta, -1.5, 1.5);
  sprintf(htitle, "cdenep%s", title);
  cdenp = new TH2D(htitle, htitle, binphi, phiL, phiU, bineta, -1.5, 1.5);

  cnump->Sumw2();
  cdenp->Sumw2();
}


void ExpCFEP::AddRealPair(double deta, double dphi, double weight)
{
  if (isoff) return;

  cnump->Fill(dphi, deta, weight);
}

void ExpCFEP::AddMixedPair(double deta, double dphi, double weight)
{
  if (isoff) return;

  cdenp->Fill(dphi, deta, weight);
}

void ExpCFEP::Write()
{
  if (isoff) return;

  cnump->Write();
  cdenp->Write();
}

