#include "ExpCF3D.h"
#include <math.h>

ExpCF3D::ExpCF3D(const char *title, 
	int binx, double minx, double maxx, 
	int biny, double miny, double maxy,
	int binz, double minz, double maxz)
{
  MakeHistos(title, binx, minx, maxx, biny, miny, maxy, binz, minz, maxz);
  isoff = false;
  useabs = false;
}

ExpCF3D::ExpCF3D(const char *title, int bin, double min, double max)
{
  MakeHistos(title, bin, min, max, bin, min, max, bin, min, max);
  isoff = false;
  useabs = false;
}

ExpCF3D::~ExpCF3D()
{ 
  if (cnum) delete cnum;
  if (cden) delete cden;
}

void ExpCF3D::SetOff(Bool_t dooff)
{
  isoff = dooff;
}

void ExpCF3D::SetUseAbs(Bool_t doabs){
  useabs = doabs;
}

void ExpCF3D::MakeHistos(const char *title, 
			 int binx, double minx, double maxx, 
			 int biny, double miny, double maxy,
			 int binz, double minz, double maxz)
{
  char htitle[1000];
  sprintf(htitle, "cnum%s", title);
  cnum = new TH3D(htitle, htitle, binx, minx, maxx, biny, miny, maxy, binz, minz, maxz);
  sprintf(htitle, "cden%s", title);
  cden = new TH3D(htitle, htitle, binx, minx, maxx, biny, miny, maxy, binz, minz, maxz);

  cnum->Sumw2();
  cden->Sumw2();
}

void ExpCF3D::AddRealPair(double qx, double qy, double qz, double weight)
{
  if (isoff) return;
  if (useabs)
    cnum->Fill(qx, qy, qz, weight);
  else
    cnum->Fill(fabs(qx), fabs(qy), fabs(qz), weight);
    
}
void ExpCF3D::AddMixedPair(double qx, double qy, double qz, double weight)
{
  if (isoff) return;
  if (useabs) 
    cden->Fill(qx, qy, qz, weight);
  else
    cden->Fill(fabs(qx), fabs(qy), fabs(qz), weight);
}

void ExpCF3D::Write()
{
  if (isoff) return;
  cnum->Write();
  cden->Write();
}
