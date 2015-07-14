#include <sys/stat.h>
#include <math.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TMath.h>

#include "CorrFctnDirectYlm.h"
#include "TPIGlobal.h"
#include "PairWeight.h"
#include "ExpCF3D.h"
#include "ExpCF1D.h"
#include "ExpCFSH.h"
#include "SourceMonitor.h"

#include <TRandom2.h>

using namespace std;

extern double ktmin, ktmax;
extern double ptmin1, ptmax1, ptmin2, ptmax2;

extern PARTICLE particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf;

extern double btmin, btmax;
extern int    pairtype;
extern int    partpid;
extern int    partpid2;
extern int    doresolution;

extern double mRL, mRS, mRO, mDT;
extern double mKO, mKS, mKL, mDE;
extern double mKOm, mKSm, mKLm, mDEm;

extern double mKStarLong, mKStarOut, mKStarSide, mKStarSigned, mRLong, mDTime, mRTrans, mROut, mRSide, mRSidePairCMS, mRLongPairCMS, mDTimePairLCMS, mROutPairCMS, mDTimePairCMS, mRStar, mBetat;
extern double mKStarLongm, mKStarOutm, mKStarSidem, mKStarSignedm, mBetatm;
extern double mKT;
extern double mKTm;
extern int twospin;
extern double pionac;

PARTICLE buf;
PARTICLE **evbuf;
PARTICLE **evbuf2;

Int_t *evtscount;
Int_t *evtscount2;


TRandom2 *mRand;


int main(int argc, char **argv)
{
  Int_t evtomix = 20;
  Double_t tcut = TMAX;
  Int_t weakpid[16] = {3334, -3334, 3312, -3312, 3322, -3322,   3112,  -3112,   3212,  -3212,   3222,  -3222,  3122, -3122,  311, -311};
  Int_t isweak;

  if (argc > 5)
    evtomix = atoi(argv[5]);

  if (argc > 6)
    tcut = atof(argv[6]);

  partpid2 = 0;

  int nbin;
  int onlyprim = 0;
  int docoulomb = 0;
  int dotrue = 0;
  int doylm = 0;
  int dosourcemon = 0;
  doresolution = 1;

  if (argc > 1)
    nbin = atoi(argv[1]);
  else {
    cout << "Usage: twointeg <kt bin> <only primaries> <do coulomb> <pair type> <events to mix> <maximum time> <do true> <do ylm> <monitor> <do resolution>" << endl;
    cout << "Pair types: " << endl << "0 - pion-pion" << endl << "1 - kaon-kaon" << endl << "2 - pion-kaon" << endl << "3 - pion-kaon unlike-sign" << endl << "4 - pion-proton" << endl << "5 - pion-proton unlike-sign" << endl << "6 - proton - Sigma+" << endl << "7 - anti-proton - Sigma+" << endl << "8 - proton - Lambda (strong)" << endl << "9 - anti-proton - Lambda (strong)" << endl << "10 - Sigma+ - Sigma+" << endl << "11 - Sigma+ - anti-Sigma+" << endl << "12 - proton-proton" << endl<< "13 - kaon+ - proton" << endl;
    exit(0);
  }

  if (argc > 2) 
    onlyprim = atoi(argv[2]);

  if (argc > 3) 
    docoulomb = atoi(argv[3]);

  if (argc > 4)
    pairtype = atoi(argv[4]);

  if (argc > 7)
    dotrue = atoi(argv[7]);

  if (argc > 8)
    doylm = atoi(argv[8]);

  if (argc > 9)
    dosourcemon = atoi(argv[9]);

  if (argc > 10)
    doresolution = atoi(argv[10]);

  btmin = 0.0;
  btmax = 1.0;

  switch(nbin) {
  case 0:
    ktmin = 0.15;
    ktmax = 0.25;
    break;
  case 1:
    ktmin = 0.25;
    ktmax = 0.35;
    break;
  case 2:
    ktmin = 0.35;
    ktmax = 0.45;
    break;
  case 3:
    ktmin = 0.45;
    ktmax = 0.6;
    break;
  case 10:
    ktmin = 0.2;
    ktmax = 0.3;
    break;
  case 11:
    ktmin = 0.3;
    ktmax = 0.36;
    break;
  case 12:
    ktmin = 0.36;
    ktmax = 0.42;
    break;
  case 13:
    ktmin = 0.42;
    ktmax = 0.48;
    break;
  case 14:
    ktmin = 0.48;
    ktmax = 0.54;
    break;
  case 15:
    ktmin = 0.54;
    ktmax = 0.60;
    break;
  case 16:
    ktmin = 0.60;
    ktmax = 0.75;
    break;
  case 17:
    ktmin = 0.75;
    ktmax = 1.0;
    break;
  case 18:
    ktmin = 1.0;
    ktmax = 1.2;
    break;
  case 19:
    ktmin = 0.2;
    ktmax = 2.0;
    break;
  case 20:
    ktmin = 0.15;
    ktmax = 0.35;
    break;
  case 21:
    ktmin = 0.35;
    ktmax = 0.55;
    break;
  case 22:
    ktmin = 0.55;
    ktmax = 0.80;
    break;
  case 30:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.3;
    btmax = 0.95;
    break;
  case 31:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.35;
    btmax = 0.5;
    break;
  case 32:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.5;
    btmax = 0.65;
    break;
  case 33:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.8;
    break;
  case 34:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.8;
    btmax = 0.95;
    break;
  case 61:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.5;
    btmax = 0.9;
    break;
  case 71:
    ktmin = 0.5;
    ktmax = 2.5;
    btmin = 0.0;
    btmax = 1.0;
    break;
  case 81:
    ktmin = 0.1;
    ktmax = 2.5;
    btmin = 0.0;
    btmax = 1.0;
    break;
  case 91:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.9;
    break;
  case 92:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.3;
    btmax = 0.99;
    break;
  case 93:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.3;
    btmax = 0.99;
    break;
  case 101:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.9;
    break;
  }

  TDatime dat;
  mRand = new TRandom2();
  mRand->SetSeed(dat.GetTime());

  switch (pairtype) {
  case 0:
    ptmin1 = PTMIN;
    ptmin2 = PTMIN;
    ptmax1 = PTMAX;
    ptmax2 = PTMAX;
    break;
  case 1:
    ptmin1 = 0.3;
    ptmin2 = 0.3;
    ptmax1 = 0.9;
    ptmax2 = 0.9;
    break;
  case 2:
  case 3:
    ptmin1 = 0.05;
    ptmin2 = 0.05;
    ptmax1 = PTMAX;
    ptmax2 = 1.5;
    break;
  case 4:
  case 5:
    ptmin1 = 0.05;
    ptmin2 = 0.3;
    ptmax1 = PTMAX;
    ptmax2 = 2.0;
    break;
  case 6:
  case 7:
    ptmin1 = 0.3;
    ptmin2 = 0.3;
    ptmax1 = 1.2;
    ptmax2 = 1.7;
    break;
  case 8:
  case 9:
    ptmin1 = 0.3;
    ptmin2 = 0.3;
    ptmax1 = 1.2;
    ptmax2 = 1.7;
    break;
  case 10:
  case 11:
    ptmin1 = 0.3;
    ptmin2 = 0.3;
    ptmax1 = 1.8;
    ptmax2 = 1.8;
    break;
  case 12:
    ptmin1 = 0.4;
    ptmin2 = 0.4;
    ptmax1 = 1.8;
    ptmax2 = 1.8;
    break;
  case 13:
    ptmin1 = 0.05;
    ptmin2 = 0.05;
    ptmax1 = 1.2;
    ptmax2 = 1.5;
    break;
  }

  if (((pairtype == 2) || (pairtype == 3)) && (nbin ==91))
    {
      ptmin1 = 0.11;
      ptmax1 = 0.5;
      ptmin2 = 0.2;
      ptmax2 = 0.9;
    }

  if (((pairtype == 4) || (pairtype == 5)) && (nbin ==92))
    {
      ptmin1 = 0.11;
      ptmax1 = 0.5;
      ptmin2 = 0.5;
      ptmax2 = 1.2;
    }

  if (((pairtype == 13) || (pairtype == 13)) && (nbin ==93))
    {
      ptmin1 = 0.3;
      ptmax1 = 0.9;
      ptmin2 = 0.5;
      ptmax2 = 1.2;
    }

  if (((pairtype == 2) || (pairtype == 3)) && (nbin ==101))
    {
      ptmin1 = 0.1;
      ptmax1 = 0.5;
      ptmin2 = 0.2;
      ptmax2 = 0.9;
//       TFile *infile1 = new TFile("piptrat.root");
//       piptrat = (TH1D *) infile1->Get("piptrat");
//       TFile *infile2 = new TFile("kptrat.root");
//       kptrat = (TH1D *) infile2->Get("kptrat");
    }

  InitializeGamow();
  if (!docoulomb)
    pionac = 0;
//   if (docoulomb) {
//     InitializeGamow();
//   }
//   else {
//     twospin = 0;
//     partpid2 = 0;
//     if (pairtype == 0) 
//       partpid = PIPID;
//     else if (pairtype == 1) 
//       partpid = KPID;
//     else if (pairtype == 10) {
//       partpid = SPPID;
//       twospin = 1;
//     }
//     else {
//       cout << "Unknown pair type " << pairtype << endl;
//       exit(0);
//     }
//   }

  TPISourceMonitorInit();
  if (!dosourcemon)
    TPISourceMonitorSetOff(kTRUE);

  evbuf = (PARTICLE **) malloc(sizeof(PARTICLE *) * evtomix);
  evtscount = (Int_t *) malloc(sizeof(Int_t) * evtomix);
  if (partpid2 != 0) {
    evbuf2 = (PARTICLE **) malloc(sizeof(PARTICLE *) * evtomix);
    evtscount2 = (Int_t *) malloc(sizeof(Int_t) * evtomix);
  }
  else {
    evbuf2 = evbuf;
    evtscount2 = evtscount;
  }

  for (int imix=0; imix<evtomix; imix++) {
    evbuf[imix] = (PARTICLE *) malloc(sizeof(PARTICLE) * 3000);
    evtscount[imix] = 0;
    if (partpid2 != 0) {
      evbuf2[imix] = (PARTICLE *) malloc(sizeof(PARTICLE) * 3000);
      evtscount2[imix] = 0;
    }
  }
  
// //   TH3D *cnumac;
// //   TH3D *cnumas;
// //   TH3D *cnuma;
// //   TH3D *cdena;

// //   TH3D *cpnumac;
// //   TH3D *cpnumas;
// //   TH3D *cpnuma;
// //   TH3D *cpdena;

// //   // 3d correlation function for non-identical particles
// //   TH3D *cnumanonid;
// //   TH3D *cdenanonid;
// //   TH3D *cnumanonidtrue;

//   TH1D *num1d;
//   TH1D *num1dqsc;
//   TH1D *num1dc;
//   TH1D *den1d;

//   TH1D *num1dl;
//   TH1D *num1dqscl;
//   TH1D *num1dcl;

//   TH1D *num1dp;
//   TH1D *den1dp;

//   TH1D *num1dn;
//   TH1D *den1dn;

//   TH1D *num1dptrue;
//   TH1D *den1dptrue;

//   TH1D *num1dntrue;
//   TH1D *den1dntrue;

  TH1D *hbetat;
  TH1D *hkt;

  ExpCF3D *cf3da;
  ExpCF3D *cf3dac;
  ExpCF3D *cf3dp;
  ExpCF3D *cf3dpc;

  if ((pairtype == 0) || (pairtype == 1) || (pairtype == 10) || (pairtype == 12)) {
    // 3D cf for identical particle - all - pure quantum weight 
    cf3da  = new ExpCF3D("a", 81, -0.2025, 0.2025);
    // 3D cf for identical particles - all - full weight
    cf3dac = new ExpCF3D("ac", 81, -0.2025, 0.2025);
    
    // 3D cf for identical particle - primaries - pure quantum weight 
    cf3dp  = new ExpCF3D("p", 81, -0.2025, 0.2025);
    // 3D cf for identical particles - primary - full weight
    cf3dpc = new ExpCF3D("pc", 81, -0.2025, 0.2025);
  }
  else {
    // 3D cf for non-identical particle - all - pure quantum weight 
    cf3da  = new ExpCF3D("a", 60, -0.2, 0.2);
    // 3D cf for non-identical particles - all - full weight
    cf3dac = new ExpCF3D("ac", 60, -0.2, 0.2);
    
    // 3D cf for non-identical particle - primaries - pure quantum weight 
    cf3dp  = new ExpCF3D("p", 60, -0.2, 0.2);
    // 3D cf for non-identical particles - primary - full weight
    cf3dpc = new ExpCF3D("pc", 60, -0.2, 0.2);
  }

// //   cnuma = new TH3D("cnuma","cnuma", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);
// //   cdena = new TH3D("cdena","cdena", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

// //   cnumac = new TH3D("cnumac","cnumac", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);
// //   cnumas = new TH3D("cnumas","cnumas", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

// //   cpnuma = new TH3D("cpnuma","cpnuma", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);
// //   cpdena = new TH3D("cpdena","cpdena", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

// //   cpnumac = new TH3D("cpnumac","cpnumac", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);
// //   cpnumas = new TH3D("cpnumas","cpnumas", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

// //   cnumanonid = new TH3D("cnumanonid","cnumanonid", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);
// //   cdenanonid = new TH3D("cdenanonid","cdenanonid", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);
// //   cnumanonidtrue = new TH3D("cnumanonidtrue","cnumanonidtrue", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);

// //   num1d = new TH1D("num1d", "num1d;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
// //   num1dqsc = new TH1D("num1dqsc", "num1dqsc;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
// //   num1dc = new TH1D("num1dc", "num1dc;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
// //   den1d = new TH1D("den1d", "den1d;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);

// //   num1dl = new TH1D("num1dl", "num1dl;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
// //   num1dqscl = new TH1D("num1dqscl", "num1dqscl;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
// //   num1dcl = new TH1D("num1dcl", "num1dcl;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);

  ExpCF1D *cf1da = new ExpCF1D("1da",200,0.0,0.4);
  ExpCF1D *cf1dp = new ExpCF1D("1dp",200,0.0,0.4);
  ExpCF1D *cf1datrue = new ExpCF1D("1datrue",200,0.0,0.4);
  
// //   num1dp = new TH1D("num1dp", "num1dp", 100, 0.0, 0.2);
// //   den1dp = new TH1D("den1dp", "den1dp", 100, 0.0, 0.2);

// //   num1dn = new TH1D("num1dn", "num1dn", 100, 0.0, 0.2);
// //   den1dn = new TH1D("den1dn", "den1dn", 100, 0.0, 0.2);

// //   num1dptrue = new TH1D("num1dptrue", "num1dptrue", 100, 0.0, 0.2);
// //   den1dptrue = new TH1D("den1dptrue", "den1dptrue", 100, 0.0, 0.2);

// //   num1dntrue = new TH1D("num1dntrue", "num1dntrue", 100, 0.0, 0.2);
// //   den1dntrue = new TH1D("den1dntrue", "den1dntrue", 100, 0.0, 0.2);

  hbetat = new TH1D("hbetat", "hbetat", 100, 0.0, 1.0);
  hkt    = new TH1D("hkt",    "hkt"   , 100, 0.0, 1.2);

// //   num1d->Sumw2();
// //   num1dqsc->Sumw2();
// //   num1dc->Sumw2();
// //   den1d->Sumw2();

// //   num1dl->Sumw2();
// //   num1dqscl->Sumw2();
// //   num1dcl->Sumw2();

// //   num1dp->Sumw2();
// //   den1dp->Sumw2();

// //   num1dn->Sumw2();
// //   den1dn->Sumw2();

// //   num1dptrue->Sumw2();
// //   den1dptrue->Sumw2();

// //   num1dntrue->Sumw2();
// //   den1dntrue->Sumw2();

// //   cnuma->Sumw2();
// //   cdena->Sumw2();

// //   cnumac->Sumw2();
// //   cnumas->Sumw2();

// //   cpnuma->Sumw2();
// //   cpdena->Sumw2();

// //   cpnumac->Sumw2();
// //   cpnumas->Sumw2();

// //   cnumanonid->Sumw2();
// //   cdenanonid->Sumw2();
// //   cnumanonidtrue->Sumw2();

  // Spherical harmonics correlation function
//   //  CorrFctnDirectYlm *cylm = new CorrFctnDirectYlm("NonIdCYlm", 3, 60, 0.0, 0.3);
//   //  CorrFctnDirectYlm *cylmtrue = new CorrFctnDirectYlm("NonIdCYlmTrue", 3, 60, 0.0, 0.3);
  ExpCFSH *cylm = new ExpCFSH("NonIdCYlm", 60, 0.0, 0.3);
  ExpCFSH *cylmtrue = new ExpCFSH("NonIdCYlmTrue", 60, 0.0, 0.3);

  // Turn the correlation functions on/off 
  // according to the switches
  if (onlyprim) {
    cf3da->SetOff(kTRUE);
    cf3dac->SetOff(kTRUE);
    cf1da->SetOff(kTRUE);
  }

  if (!docoulomb) {
    cf3dac->SetOff(kTRUE);
    cf3dpc->SetOff(kTRUE);
  }

  if (!dotrue) {
    cf1datrue->SetOff(kTRUE);
    cylmtrue->SetOff(kTRUE);
  }

  if (!doylm) {
    cylm->SetOff(kTRUE);
    cylmtrue->SetOff(kTRUE);
  }

  TChain *chn = new TChain("particles");
  int pcount = 0;
  
  char fname[100];
  int nfile = 1;
  struct stat fbuf;

  sprintf(fname, "event1.root");
  while (!stat(fname, &fbuf)) {
    cout << "Adding " << fname << endl;
    chn->Add(fname);
    nfile++;
    sprintf(fname, "event%i.root", nfile);
  }
  
  Int_t npart = (int) chn->GetEntries();

  //it must be initialized
  chn->SetBranchAddress("part",&buf);
  Int_t curev = -1;
  Int_t eviter = -1;
  Double_t pt, rap, peta;
  double coulombweight;
  double quantumweight;
  double iscorrelated;

  for (int iter=0; iter<npart; iter++) {
    chn->GetEntry(iter);
    if (((buf.pid == partpid) || ((buf.pid == partpid2) && (partpid2 !=0))) && (buf.t<tcut)) {
      pt = hypot(buf.px, buf.py);
      if ((((pt>ptmin1) && (pt<ptmax1) && (buf.pid == partpid)) ||
	   ((pt>ptmin2) && (pt<ptmax2) && (buf.pid == partpid2) && (partpid2 != 0))) &&
	  ((onlyprim == 0) || ((buf.fatherpid == partpid) || ((buf.fatherpid == partpid2) && (partpid2 != 0))))) {

	rap = 0.5*log((buf.e+buf.pz)/(buf.e-buf.pz));
	peta = -TMath::Log(TMath::Tan(TMath::ATan2(pt, buf.pz)/2.0));
	
	// Check if it is coming from a weak decay
	isweak = 0;
	for (int witer=0; witer<16; witer++) 
	  if (buf.fatherpid == weakpid[witer]) isweak = 1;
	      
	if ((fabs(rap) < ABSRAP) &&
	    ((fabs(peta) < ETAABS) || (nbin < 10) || (nbin>19)) &&
	    ((!isweak) || (onlyprim == 2)))
	  {
	    if (buf.pid == partpid) {
	      evbuf[eviter][evtscount[eviter]].x = buf.x;
	      evbuf[eviter][evtscount[eviter]].y = buf.y;
	      evbuf[eviter][evtscount[eviter]].z = buf.z;
	      evbuf[eviter][evtscount[eviter]].t = buf.t;
	      evbuf[eviter][evtscount[eviter]].px = buf.px;
	      evbuf[eviter][evtscount[eviter]].py = buf.py;
	      evbuf[eviter][evtscount[eviter]].pz = buf.pz;
	      evbuf[eviter][evtscount[eviter]].e = buf.e;
	      evbuf[eviter][evtscount[eviter]].pid = buf.pid;
	      evbuf[eviter][evtscount[eviter]].fatherpid = buf.fatherpid;
	      evbuf[eviter][evtscount[eviter]].rootpid = buf.rootpid;
	      evbuf[eviter][evtscount[eviter]].mass = buf.mass;
	      evbuf[eviter][evtscount[eviter]].event = buf.event;
	      evtscount[eviter]++;
	    }
	    else if ((buf.pid == partpid2) && (partpid2 != 0)) {
	      evbuf2[eviter][evtscount2[eviter]].x = buf.x;
	      evbuf2[eviter][evtscount2[eviter]].y = buf.y;
	      evbuf2[eviter][evtscount2[eviter]].z = buf.z;
	      evbuf2[eviter][evtscount2[eviter]].t = buf.t;
	      evbuf2[eviter][evtscount2[eviter]].px = buf.px;
	      evbuf2[eviter][evtscount2[eviter]].py = buf.py;
	      evbuf2[eviter][evtscount2[eviter]].pz = buf.pz;
	      evbuf2[eviter][evtscount2[eviter]].e = buf.e;
	      evbuf2[eviter][evtscount2[eviter]].pid = buf.pid;
	      evbuf2[eviter][evtscount2[eviter]].fatherpid = buf.fatherpid;
	      evbuf2[eviter][evtscount2[eviter]].rootpid = buf.rootpid;
	      evbuf2[eviter][evtscount2[eviter]].mass = buf.mass;
	      evbuf2[eviter][evtscount2[eviter]].event = buf.event;
	      evtscount2[eviter]++;
	    }
	  }
      }
    }

    if (buf.event != curev) {
      if (eviter>-1) {
	// Mix particles 
	for (int mixiter=0; mixiter<evtomix; mixiter++) {
	  for (int fiter=0; fiter<evtscount[eviter]; fiter++) {
	    for (int siter=0; siter<evtscount2[mixiter]; siter++) {
	      if ((partpid2 == 0) && (mixiter==eviter) &&  (siter>=fiter)) break; 
	      PairKinematics(evbuf[eviter][fiter], 
			     evbuf2[mixiter][siter]);
	      if ((mKTm>ktmin) && (mKTm<ktmax)) {
		cf3da->AddMixedPair(mKOm, mKSm, mKLm, 1.0);
		if (partpid2 == 0)
		  cf1da->AddMixedPair(fabs(mKStarSignedm)*2.0, 1.0);
		else 
		  cf1da->AddMixedPair(fabs(mKStarSignedm), 1.0, mKStarOutm>0);
		cylm->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		
		if (partpid2 != 0) {
		  hbetat->Fill(mBetatm);
		  hkt->Fill(mKTm);
		}

		if (docoulomb) {
		  cf3dac->AddMixedPair(mKOm, mKSm, mKLm, 1.0);

		  quantumweight = GetQuantum();
		  cf3da->AddRealPair(mKOm, mKSm, mKLm, quantumweight);

		  coulombweight = GetFull();
		    
		  if (partpid2 == 0)
		    cf1da->AddRealPair(fabs(mKStarSignedm)*2.0, coulombweight);
		  else
		    cf1da->AddRealPair(fabs(mKStarSignedm), coulombweight, mKStarOutm>0);
		  
		  cf3dac->AddRealPair(mKOm, mKSm, mKLm, coulombweight);
		  cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
		  
		  if (mixiter == eviter) {
		    cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
		    if (partpid2 == 0)
		      cf1datrue->AddRealPair(fabs(mKStarSignedm)*2.0, coulombweight);
		    else
		      cf1datrue->AddRealPair(fabs(mKStarSignedm), coulombweight, mKStarOutm>0);
		  }
		  else {
		    cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		    if (partpid2 == 0)
		      cf1datrue->AddMixedPair(fabs(mKStarSignedm)*2.0, 1.0);
		    else
		      cf1datrue->AddMixedPair(fabs(mKStarSignedm), 1.0, mKStarOutm>0);
		  }
		}
		else {
		  quantumweight = GetQuantum();
		  
		  cf1da->AddRealPair(fabs(mKStarSignedm)*2.0, quantumweight);
		  cf3da->AddRealPair(mKOm, mKSm, mKLm, quantumweight);
		  cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, quantumweight);

		  if (mixiter == eviter) {
		    cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, quantumweight);
		    if (partpid2 == 0)
		      cf1datrue->AddRealPair(fabs(mKStarSignedm)*2.0, quantumweight);
		    else
		      cf1datrue->AddRealPair(fabs(mKStarSignedm), quantumweight, mKStarOutm>0);
		  }
		  else {
		    cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		    if (partpid2 == 0)
		      cf1datrue->AddMixedPair(fabs(mKStarSignedm)*2.0, 1.0);
		    else
		      cf1datrue->AddMixedPair(fabs(mKStarSignedm), 1.0, mKStarOutm>0);
		  }
		}
	      }
	    }
	  }
	  if (partpid2 == 0) continue;
	  for (int fiter=0; fiter<evtscount[mixiter]; fiter++) {
	    for (int siter=0; siter<evtscount2[eviter]; siter++) {
	      PairKinematics(evbuf[mixiter][fiter], 
			     evbuf2[eviter][siter]);
	      if ((mKTm>ktmin) && (mKTm<ktmax)) {
		cf3dac->AddMixedPair(mKOm, mKSm, mKLm, 1.0);
		cf1da->AddMixedPair(fabs(mKStarSignedm), 1.0, mKStarOutm>0);
		cylm->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		
		if (partpid2 != 0) {
		  hbetat->Fill(mBetatm);
		  hkt->Fill(mKTm);
		}
		
		if (docoulomb) {
		  coulombweight = GetFull();
		  
		  cf1da->AddRealPair(fabs(mKStarSignedm), coulombweight, mKStarOutm>0);
		  
		  cf3dac->AddRealPair(mKOm, mKSm, mKLm, coulombweight);
		  cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
		  
		  if (mixiter == eviter) {
		    cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
		    cf1datrue->AddRealPair(fabs(mKStarSignedm), coulombweight, mKStarOutm>0);
		  }
		  else {
		    cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
		    cf1datrue->AddMixedPair(fabs(mKStarSignedm), 1.0, mKStarOutm>0);
		  }
		}
	      }
	    }
	  }
	}
      }
        

// // 	if (partpid2 == 0) {
// // 	  // Mix different-event particles;
// // 	  // for identical particle correlations
// // 		  if ((mKTm>ktmin) && (mKTm<ktmax)) {
// // 		    cf3da->AddMixedPair(mKOm, mKSm, mKLm, 1.0);
// // 		    //		    cpdena->Fill(mKStarOutm, mKStarSidem, mKStarLongm, 
// // 		    //				1.0);
// // 		    cf1da->AddMixedPair(fabs(mKStarSignedm)*2.0, 1.0);
// // 		    cylm->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);
// // 		    cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1.0);

// // 		    if (docoulomb)
// // 		      //		      if ((fabs(mKO)<0.4) && (fabs(mKS)<0.4) && (fabs(mKL)<0.4)) 
// // 		      {


// // 			cnumas->Fill(mKOm, mKSm, mKLm, coulombweight);
// // 			cpnumas->Fill(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// // 			num1dqsc->Fill(fabs(mKStarSignedm)*2.0, coulombweight);
// // 			num1dc->Fill(fabs(mKStarSignedm)*2.0, coulqscpart);
// //  			cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// //  			cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// // 		      }
// // 		  }
		  
// // 		}
// // 	      }
// // 	    }
// // 	  }
// // 	}
	
// 		  // Mix different-event particles
// 		  // for non-identical particle correlations
// // 		  else {
// // 	  for (int fiter=0; fiter<evtscount[eviter]; fiter++) {
// // 	    for (int mixiter=0; mixiter<evtomix; mixiter++) {
// // 	      //	      if (mixiter != eviter) {
// // 	      for (int siter=0; siter<evtscount2[mixiter]; siter++) {
// // 		PairKinematics(evbuf[eviter][fiter], 
// // 			       evbuf2[mixiter][siter]);
// // 		if ((mKTm>ktmin) && (mKTm<ktmax) && (mBetatm>btmin) && (mBetatm<btmax)) {
// // 		  cdena->Fill(mKOm, mKSm, mKLm, 
// // 			      1.0);
// // 		  den1d->Fill(fabs(mKStarSignedm)*2.0, 1.0);


// // 		  if (mKStarOutm > 0.0) 
// // 		    den1dp->Fill(fabs(mKStarSignedm), 1.0);
// // 		  else
// // 		    den1dn->Fill(fabs(mKStarSignedm), 1.0);

// // 		  cdenanonid->Fill(mKStarOutm, mKStarSidem, mKStarLongm, 1);
		  
// //  		  cylm->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1);
// //  		  cylmtrue->AddMixedPair(mKStarOutm, mKStarSidem, mKStarLongm, 1);

// // 		  if (mKStarOutm > 0.0) {
// //  		    if (mixiter != eviter) 
// // 		      den1dptrue->Fill(fabs(mKStarSignedm), 1.0);
// // 		  }
// // 		  else {
// //  		    if (mixiter != eviter) 
// // 		      den1dntrue->Fill(fabs(mKStarSignedm), 1.0);
// // 		  }
// // 		  if (docoulomb) {
// // 		    //		    if ((fabs(mKO)<0.2) && (fabs(mKS)<0.2) && (fabs(mKL)<0.2)) {
// // 		    if (fabs(pionac) > 0.1) {
// // 		      if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
// // 			coulombweight = 1.0;
// // 		      else
// // 			coulombweight = GetCoulomb();
// // 		      num1d->Fill(fabs(mKStarSignedm), coulombweight); 
// // 		      cnumanonid->Fill(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);

// //  		      cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
		      
// // 		      if (mKStarOutm > 0.0)
// // 			num1dp->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      else
// // 			num1dn->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      if (mixiter == eviter) {
// // 			cnumanonidtrue->Fill(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// //  			cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// // 			if (mKStarOutm > 0.0)
// // 			  num1dptrue->Fill(fabs(mKStarSignedm), coulombweight);
// // 			else
// // 			  num1dntrue->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      }
// // 		    }
// // 		    else {
// // 		      coulombweight = GetStrong();
// // 		      num1d->Fill(fabs(mKStarSignedm)*2.0, coulombweight); 
// // 		      cnumanonid->Fill(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);

// //  		      cylm->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);

// // 		      if (mKStarOutm > 0.0)
// // 			num1dp->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      else
// // 			num1dn->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      if (mixiter == eviter) {
// // 			cnumanonidtrue->Fill(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// //  			cylmtrue->AddRealPair(mKStarOutm, mKStarSidem, mKStarLongm, coulombweight);
// // 			if (mKStarOutm > 0.0)
// // 			  num1dptrue->Fill(fabs(mKStarSignedm), coulombweight);
// // 			else
// // 			  num1dntrue->Fill(fabs(mKStarSignedm), coulombweight);
// // 		      }
// // 		    }
// // 		  }
// // 		}
// // 		//		}
		
// // 	      }
// // 	    }
// // 	  }
// // 	}
// //       }
      
      curev = buf.event;
      eviter = (eviter + 1) % evtomix;
      
      cout << "Event: "<< curev << " particle 1 count " << evtscount[eviter] <<  endl;
      evtscount[eviter] = 0;
      if (partpid2 != 0) {
	cout << "Event: "<< curev << " particle 2 count " << evtscount2[eviter] <<  endl;
	evtscount2[eviter] = 0;
      }      
    }
  }
  
  char bufs[100];
  if (doresolution) {
    if (pairtype == 0)
      sprintf(bufs, "outfilecf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 1)
      sprintf(bufs, "outfilekkcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 2)
      sprintf(bufs, "outfilepikcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 3)
      sprintf(bufs, "outfilepikulcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 4)
      sprintf(bufs, "outfilepipcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 5)
      sprintf(bufs, "outfilepipulcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 6)
      sprintf(bufs, "outfilepspulcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 7)
      sprintf(bufs, "outfileapspulcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 8)
      sprintf(bufs, "outfileplacf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 9)
      sprintf(bufs, "outfileaplacf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 10)
      sprintf(bufs, "outfilespspcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 11)
      sprintf(bufs, "outfilespspulcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 12)
      sprintf(bufs, "outfileppcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 13)
      sprintf(bufs, "outfilekpcf%i%s_mr.root", nbin, onlyprim ? "p" : "a");
  }
  else {
    if (pairtype == 0)
      sprintf(bufs, "outfilecf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 1)
      sprintf(bufs, "outfilekkcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 2)
      sprintf(bufs, "outfilepikcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 3)
      sprintf(bufs, "outfilepikulcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 4)
      sprintf(bufs, "outfilepipcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 5)
      sprintf(bufs, "outfilepipulcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 6)
      sprintf(bufs, "outfilepspulcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 7)
      sprintf(bufs, "outfileapspulcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 8)
      sprintf(bufs, "outfileplacf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 9)
      sprintf(bufs, "outfileaplacf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 10)
      sprintf(bufs, "outfilespspcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 11)
      sprintf(bufs, "outfilespspulcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 12)
      sprintf(bufs, "outfileppcf%i%s.root", nbin, onlyprim ? "p" : "a");
    else if (pairtype == 13)
      sprintf(bufs, "outfilekpcf%i%s.root", nbin, onlyprim ? "p" : "a");
  }

  TFile *ofile = new TFile(bufs, "RECREATE");
  ofile->cd();
  
  cylm->Write();
  cylmtrue->Write();

  cf1da->Write();
  cf1datrue->Write();

  cf3da->Write();
  cf3dac->Write();
  
  cf3dp->Write();
  cf3dpc->Write();
  
  hbetat->Write();
  hkt->Write();

}
