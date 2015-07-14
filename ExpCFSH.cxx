#include "ExpCFSH.h"

ExpCFSH::ExpCFSH(const char *title, int bin, double min, double max)
{
  isoff = false;
  MakeHistos(title, bin, min, max);
}

ExpCFSH::~ExpCFSH()
{
  if (cfun) delete cfun;
}
  
void ExpCFSH::SetOff(Bool_t dooff)
{
  isoff = dooff;
}

void ExpCFSH::MakeHistos(const char *title, int bin, double min, double max)
{
  cfun = new CorrFctnDirectYlm(title, 3, bin, min, max);
}


void ExpCFSH::AddRealPair(double qx, double qy, double qz, double weight)
{
  if (isoff) return;

  cfun->AddRealPair(qx, qy, qz, weight);
}

void ExpCFSH::AddMixedPair(double qx, double qy, double qz, double weight)
{
  if (isoff) return;

  cfun->AddMixedPair(qx, qy, qz, weight);
}

void ExpCFSH::Write()
{
  if (isoff) return;

  cfun->Finish();
  cfun->Write();
}

