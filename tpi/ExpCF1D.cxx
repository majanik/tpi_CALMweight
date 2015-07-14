#include "ExpCF1D.h"

ExpCF1D::ExpCF1D(const char *title, int bin, double min, double max)
{
  isoff = false;
  MakeHistos(title, bin, min, max);
}

ExpCF1D::~ExpCF1D()
{
  if (cnump) delete cnump;
  if (cdenp) delete cdenp;
  if (cnumn) delete cnumn;
  if (cdenn) delete cdenn;
}
  
void ExpCF1D::SetOff(Bool_t dooff)
{
  isoff = dooff;
}

void ExpCF1D::MakeHistos(const char *title, int bin, double min, double max)
{
  char htitle[2000];

  sprintf(htitle, "cnump%s", title);
  cnump = new TH1D(htitle, htitle, bin, min, max);
  sprintf(htitle, "cdenp%s", title);
  cdenp = new TH1D(htitle, htitle, bin, min, max);
  sprintf(htitle, "cnumn%s", title);
  cnumn = new TH1D(htitle, htitle, bin, min, max);
  sprintf(htitle, "cdenn%s", title);
  cdenn = new TH1D(htitle, htitle, bin, min, max);

  cnump->Sumw2();
  cnumn->Sumw2();
  cdenp->Sumw2();
  cdenn->Sumw2();
}


void ExpCF1D::AddRealPair(double qv, double weight, int sw)
{
  if (isoff) return;

  if (sw) 
    cnump->Fill(qv, weight);
  else
    cnumn->Fill(qv, weight);
}

void ExpCF1D::AddMixedPair(double qv, double weight, int sw)
{
  if (isoff) return;

  if (sw)
    cdenp->Fill(qv, weight);
  else
    cdenn->Fill(qv, weight);
}

void ExpCF1D::Write()
{
  if (isoff) return;

  cnumn->Write();
  cdenn->Write();

  if ((cdenp->GetEntries()>0) || (cnump->GetEntries()>0)) {
    // Create also the sum hist
    char htitle[2000];
    TH1D *cnums = new TH1D(*cnumn);
    cnums->Add(cnump);
    sprintf(htitle, "cnum");
    strcat(htitle, cnump->GetName()+5);
    cnums->SetName(htitle);
    cnums->SetTitle(htitle);

    TH1D *cdens = new TH1D(*cdenn);
    cdens->Add(cdenp);
    sprintf(htitle, "cden");
    strcat(htitle, cdenp->GetName()+5);
    cdens->SetName(htitle);
    cdens->SetTitle(htitle);
  
    cnump->Write();
    cdenp->Write();

    cnums->Write();
    cdens->Write();

  }

}

