#ifndef _TPF_ExpCF3D_
#define _TPF_ExpCF3D_

#include <TH3D.h>

class ExpCF3D
{
 public:
  ExpCF3D(const char *title, int bin, double min, double max);
  ExpCF3D(const char *title, 
	  int binx, double minx, double maxx, 
	  int biny, double miny, double maxy,
	  int binz, double minz, double maxz);
  ~ExpCF3D();
  
  void AddRealPair(double qx, double qy, double qz, double weight);
  void AddMixedPair(double qx, double qy, double qz, double weight);
  
  void Write();
  void SetOff(Bool_t dooff);
  void SetUseAbs(Bool_t doabs);

 private:

  int isoff;
  int useabs;
  void MakeHistos(const char *title, 
	  int binx, double minx, double maxx, 
	  int biny, double miny, double maxy,
	  int binz, double minz, double maxz);

  TH3D *cnum;
  TH3D *cden;
};

#endif
