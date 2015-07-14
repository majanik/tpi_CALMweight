#include <math.h>
#include <TMath.h>
#include <TRandom2.h>

#include "SourceMonitor.h"
#include "PairWeight.h"

// external variables
extern TRandom2 *mRand;

//PARTICLE particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf;
//PARTICLE particle1prfm, particle1lcmsm, particle2prfm, particle2lcmsm;
ParticleCoor particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf;
ParticleCoor particle1prfm, particle1lcmsm, particle2prfm, particle2lcmsm;

struct MomResParameters{
  double PhiA;
  double PhiB;
  double PhiAlpha;
  double PhiMu;
  double ThetaA;
  double ThetaB;
  double ThetaAlpha;
  double PMeanA;
  double PMeanB;
  double PMeanAlpha;
  double PResA;
  double PResB;
  double PResAlpha;
  double PResC;
};

double mKStarLong, mKStarOut, mKStarSide, mKStarSigned, mRLong, mDTime, mRTrans, mROut, mRSide, mRSidePairCMS, mRLongPairCMS, mDTimePairLCMS, mROutPairCMS, mDTimePairCMS, mRStar, mBetat;
double mKStarLongm, mKStarOutm, mKStarSidem, mKStarSignedm, mBetatm;
double pairphi;
double pairphim;

double mRL, mRS, mRO, mDT;
double mKO, mKS, mKL, mDE;
double mKOm, mKSm, mKLm, mDEm;
double mKOL, mKSL, MKLL;
double deta, dphi;
double phiLow, phiUp;
double ev2weight;

double mROS, mRSS, mRLS, mRSt;
double mKT;
double mKTm;
double ktmin, ktmax;
double ptmin1, ptmax1, ptmin2, ptmax2;

double btmin, btmax;
int    uselcms;
int    pairtype;
int    partpid;
int    partpid2;
int    doresolution;

double mMomResScale;
MomResParameters* mMomRes1;
MomResParameters* mMomRes2;  
static int pcount = 0;

// end external variables

struct _dcomplex 
{
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;

dcomplex d0s;
dcomplex f0s;
dcomplex d0t;
dcomplex f0t;

double gamov[2000];
long double pionac;
long double oneoveracsq;
long double twopioverac;
long double coulqscpart;     
int  twospin;
int  writegrps;
long double euler;
long double f0;
long double d0;

MomResParameters* GetPionMomResParameters(){
  MomResParameters* t = new MomResParameters;
  t->PhiA = 0.00201;
  t->PhiB = 0.001018;
  t->PhiAlpha = -1.274;
  t->ThetaA = 0.000908;
  t->ThetaB = 0.001255;
  t->ThetaAlpha = -1.141;
  t->PResA = 0.01074;
  t->PResB = 0.001918;
  t->PResAlpha = -0.9895;
  t->PResC = 0.009706;
  t->PMeanA = 0.; // no energy loss (Kalman filter)
  t->PMeanB = 0.;
  t->PMeanAlpha = 0.;
  return t;
}

MomResParameters* GetKaonMomResParameters(){
  MomResParameters* t = new MomResParameters;
  t->PResA = 0.01981;
  t->PResB = 0.001371;
  t->PResC = 0.0;
  t->PResAlpha = -2.112;
  t->PhiA = 0.001791;
  t->PhiB = 0.001319;
  t->PhiAlpha = -1.686;
  t->ThetaA = 0.0005202;
  t->ThetaB = 0.001752;
  t->ThetaAlpha = -1.352;
  t->PMeanA = -0.004136;
  t->PMeanB = 0.003511;
  t->PMeanAlpha = -1.192;
  return t;
}

MomResParameters* GetProtonMomResParameters(){
  MomResParameters* t = new MomResParameters;
  t->PhiA = 0.0006575;
  t->PhiB = 0.002813;
  t->PhiAlpha = -1.583;
  t->ThetaA = 0.0002846;
  t->ThetaB = 0.002458;
  t->ThetaAlpha = -1.475;
  t->PResA = 0.7989;
  t->PResB = -2.128;
  t->PResAlpha = 0.5282;
  t->PResC = 1.41;
  t->PMeanA = -0.006509;
  t->PMeanB = 0.008757;
  t->PMeanAlpha = -1.373;
  return t;
}

void InitializeGamow()
{
  twospin = 0;
  mMomResScale = 1.0;
  euler = 0.577215665;
  f0 = 7.77 / 0.197327;
  d0 = 2.77 / 0.197327;

  d0s.re = 0.0;
  d0s.im = 0.0;
  d0t.re = 0.0;
  d0t.im = 0.0;

  f0s.re = 0.0;
  f0s.im = 0.0;
  f0t.re = 0.0;
  f0t.im = 0.0;

  mMomRes1 = 0;
  mMomRes2 = 0;

  if (pairtype == 0) {
    pionac = 387.5 / 0.197327;
    partpid = PIPID;
    partpid2 = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  }
  else if (pairtype == 1) {
    pionac = 109.55 / 0.197327;
    partpid = KPID;
    partpid2 = 0;
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
  }
  else if (pairtype == 2) {
    pionac = 248.519 / 0.197327;
    partpid = PIPID;
    partpid2 = KPID;
    writegrps = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
  }
  else if (pairtype == 3) {
    pionac = -248.519 / 0.197327;
    partpid = -PIPID;
    partpid2 = KPID;
    writegrps = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
  }
  else if (pairtype == 4) {
    pionac = 222.564 / 0.197327;
    partpid = PIPID;
    partpid2 = PPID;
    writegrps = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
  }
  else if (pairtype == 5) {
    pionac = -222.564 / 0.197327;
    partpid = -PIPID;
    partpid2 = PPID;
    writegrps = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
  }
  else if (pairtype == 6) {
    pionac = 51.5553 / 0.197327;
    partpid = PPID;
    partpid2 = SPPID;
    twospin = 1;
    // Rijken nucl-th/9901028 model a
//     d0s.re = 3.16 / 0.197327;
//     f0s.re = 4.35 / 0.197327;
//     d0t.re = -59.48 / 0.197327;
//     f0t.re = 0.25 / 0.197327;

    // Rijken, Nagels, Yamamoto ESC08 model PTPS No. 185, 2010
    f0s.re = 3.85 / 0.197327;
    d0s.re = 3.40 / 0.197327;
    f0t.re = -0.62 / 0.197327;
    d0t.re = -2.13 / 0.197327;
  }
  else if (pairtype == 7) {
    pionac = -51.5553 / 0.197327;
    partpid = -PPID;
    partpid2 = SPPID;

    // from antiproton - lambda
    f0s.re = -1.20 /0.197327;
    f0s.im = 0.37  /0.197327;
    f0t.re = -1.20 /0.197327;
    f0t.im = 0.37  /0.197327;
  }
  else if (pairtype == 8) {
    pionac = 0.0;
    partpid = PPID;
    partpid2 = LAPID;
      
//    f0s.re = 2.88 /0.197327;
//    d0s.re = 2.92 /0.197327;
//    f0t.re = 1.66 /0.197327;
//    d0t.re = 3.78 /0.197327;
      
    // Rijken, Nagels, Yamamoto ESC08 model PTPS No. 185, 2010
    f0s.re = 2.70 / 0.197327;
    d0s.re = 2.97 / 0.197327;
    f0t.re = 1.65 / 0.197327;
    d0t.re = 3.63 / 0.197327;
      
  }
  else if (pairtype == 9) {
    pionac = 0.0;
    partpid = -PPID;
    partpid2 = LAPID;
      
    f0s.re = -1.20 /0.197327;
    f0s.im = 0.37  /0.197327;
    f0t.re = -1.20 /0.197327;
    f0t.im = 0.37  /0.197327;
  }
  else if (pairtype == 10) {
    pionac = 45.4709 / 0.197327;
    partpid = SPPID;
    partpid2 = 0;
    twospin = 1;
      
    //no data from ESC08
  }
  else if (pairtype == 11) {
    pionac = -45.4709 / 0.197327;
    partpid = SPPID;
    partpid2 = -SPPID;
      
    //no data from ESC08
  }
  else if (pairtype == 12) {
    pionac = 57.63975274 / 0.197327;
    partpid = PPID;
    partpid2 = 0;
    twospin = 1;
      
//    d0s.re = 2.77 / 0.197327;
//    f0s.re = 7.77 / 0.197327;
    d0t.re = 1.7  / 0.197327;
    f0t.re = -5.4 / 0.197327;
      
    // Rijken, Nagels, Yamamoto ESC08 model PTPS No. 185, 2010
    f0s.re = 7.771 / 0.197327;
    d0s.re = 2.754 / 0.197327;
//    f0t.re = ? / 0.197327; no data from ESC08
//    d0t.re = ? / 0.197327;
  }
  else if (pairtype == 13) {
    pionac = 83.59432006 / 0.197327;
    partpid = KPID;
    partpid2 = PPID;
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
  }
  else if (pairtype == 14) {
    pionac = -83.59432006 / 0.197327;
    partpid = -KPID;
    partpid2 = PPID;
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
  }
  else if (pairtype == 15) {
    pionac = -387.5 / 0.197327;
    partpid = PIPID;
    partpid2 = -PIPID;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  }
  else if (pairtype == 16) {
    pionac = -387.5 / 0.197327;
    partpid = PIPID;
    partpid2 = -PIPID;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  }
  else if (pairtype == 17) {
    pionac = 48.5 / 0.197327;
    partpid = LAPID;
    partpid2 = 0;

    // Rijken, Nagels, Yamamoto ESC08 model PTPS No. 185, 2010
    f0s.re = 0.88 / 0.197327;
    d0s.re = 4.34 / 0.197327;
//    f0t.re = ? / 0.197327; no data from ESC08
//    d0t.re = ? / 0.197327;
  }
  else if (pairtype == 18) {
      pionac = 49.4 / 0.197327;
      partpid = PPID;
      partpid2 = XIPID;
      twospin = 1;
     
      // Rijken, Nagels, Yamamoto ESC08 model PTPS No. 185, 2010
      f0s.re = -0.58 / 0.197327;
      d0s.re = -2.71 / 0.197327;
      f0t.re = -3.49 / 0.197327;
      d0t.re =  0.60 / 0.197327;
  }
  else if (pairtype == 19) {
    pionac =  -109.55 / 0.197327;
    partpid = KPID;
    partpid2 = -KPID;
    mMomRes1 = GetKaonMomResParameters();
    mMomRes2 = GetKaonMomResParameters();
  }
  else if (pairtype == 20) {
    pionac =  -57.63975274 / 0.197327;
    partpid = PPID;
    partpid2 = -PPID;
    mMomRes1 = GetProtonMomResParameters();
    mMomRes2 = GetProtonMomResParameters();
  }
  else if (pairtype == 21) {
    pionac = 0;//387.5 / 0.197327;
    partpid = PI0PID;
    partpid2 = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  }
  else if (pairtype == 22) {
    pionac = 387.5 / 0.197327;
    partpid = -PIPID;
    partpid2 = 0;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  }  
  else if (pairtype == 23) {
    pionac = 387.5 / 0.197327;
    partpid = PIPID;
    partpid2 = PIPID;
    mMomRes1 = GetPionMomResParameters();
    mMomRes2 = GetPionMomResParameters();
  } 
  else {
    cout << "Unknown pair type " << pairtype << endl;
    exit(0);
  }

  if (!doresolution) {
    mMomRes1 = 0;
    mMomRes2 = 0;
  }

  oneoveracsq = 1.0/(pionac * pionac);
  twopioverac = 2.0*TMath::Pi()/pionac;
  double tpaoverk;
  
  //Zajmuje duzo czasu, Jeremi mowi zeby wywalic :D
  // for (int iter=0; iter<2000; iter++) {
  //   tpaoverk = twopioverac/(iter*0.0002 + 0.0001);
  //   gamov[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
  // }
}

double Gamow(double arg)
{
  //   if (arg<0.0001)
  //     return 0.0;
  //   if (arg > 0.4)
  long double eta = twopioverac / arg;
  return (eta) * 1.0/(exp(eta) - 1.0);
  //   int bin = arg / 0.0002;
  //   double b1 = bin*0.0002 + 0.0001;
  //   double b2 = bin*0.0002 + 0.0003;
  //   return ((gamov[bin+1] * (arg - b1) + gamov[bin] * (b2 - arg)) / 0.0002);
}

dcomplex mult(dcomplex arga, dcomplex argb)
{
  dcomplex res;
  
  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}

dcomplex conj(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re;
  res.im = -arg.im;

  return res;
}

dcomplex mult(dcomplex arga, long double argb)
{
  dcomplex res;
  
  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

long double modl2(dcomplex arg)
{
  return arg.re*arg.re + arg.im*arg.im;
}

long double modl(dcomplex arg)
{
  return hypot(arg.re, arg.im);
}

dcomplex invr(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re/modl2(arg);
  res.im = - arg.im/modl2(arg);

  return res;
}

// Calculates the confluent hypergeometric function F
// from single orientation of cos(theta*)
// For non-symmetrized wave-function (non-identical particles)
void GetFFsingle(dcomplex *ffp, int sign = 1)
{
  double comprep[COULOMBSTEPS];
  double compimp[COULOMBSTEPS];
  double eta, ksip;
  dcomplex alfa, zetp;
  
  int nsteps;

  double kstar = fabs(mKStarSigned);
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  double coskr = sign * kstrst/(fabs(mKStarSigned) * mRSt);

  if (kstar*mRSt*(1+coskr) > 15.0)
    nsteps = 170;
  else if (kstar*mRSt*(1+coskr) > 10.0)
    nsteps = 45;
  else if (kstar*mRSt*(1+coskr) > 5.0)
    nsteps = 35;
  else
    nsteps = 15;

  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp;
  double tcomp;
  dcomplex sump;
  dcomplex fcmult;

  double rad = mRSt;

  ksip = kstar*rad*(1+coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);

    sump = mult(sump, 1.0/(tcomp*tcomp));
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    tcomp *= (istep+1);

    if ((sump.re*sump.re+sump.im*sump.im) < 1.0e-14) {
      nsteps = istep;
      break;
    }
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  
}

// Calculates the confluent hypergeometric function
// For two orientations of cos(theta*) 
// For symmetrized wave-function (identical particles)
void GetFFdouble(dcomplex *ffp, dcomplex *ffm)
{
  long double comprep[COULOMBSTEPS];
  long double compimp[COULOMBSTEPS];
  long double comprem[COULOMBSTEPS];
  long double compimm[COULOMBSTEPS];
  long double eta, ksip, ksim;
  dcomplex alfa, zetp, zetm;
  
  int nsteps;

  long double kstar = fabs(mKStarSigned);
  long double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double coskr = kstrst/(fabs(mKStarSigned) * mRSt);

  if ((kstar*mRSt*(1+coskr) < 5.0) &&
      (kstar*mRSt*(1-coskr) < 5.0))
    nsteps = 25;
  else if ((kstar*mRSt*(1+coskr) < 10.0) &&
	   (kstar*mRSt*(1-coskr) < 10.0))
    nsteps = 45;
  else if ((kstar*mRSt*(1+coskr) < 15.0) &&
	   (kstar*mRSt*(1-coskr) < 15.0))
    nsteps = 110;
  else
    nsteps = 150;
  
  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp, scompm;
  long double tcomp;
  dcomplex sump, summ;
  dcomplex fcmult;

  long double rad = mRSt;

  ksip = kstar*rad*(1+coskr);
  ksim = kstar*rad*(1-coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  zetm.re = 0.0;
  zetm.im = ksim;

  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  scompm.re = 1.0; 
  scompm.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);
    summ = mult(fcomp, scompm);

    sump = mult(sump, 1.0/(tcomp*tcomp));
    summ = mult(summ, 1.0/(tcomp*tcomp));
	
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
      
      comprem[istep] = summ.re;
      compimm[istep] = summ.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
      
      comprem[istep] = comprem[istep-1] + summ.re;
      compimm[istep] = compimm[istep-1] + summ.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    scompm = mult(scompm, zetm);
    tcomp *= (istep+1);
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  ffm->re = comprem[nsteps-1];
  ffm->im = compimm[nsteps-1];
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized)
// with Coulomb interaction included
double GetQuantumCoulomb()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  int ccase = 0;
  int wavesign = 1;

  if (twospin == 1) {
    if (pcount == 3)
      wavesign = 1;
    else
      wavesign = -1;
    pcount++;
    if (pcount == 4) pcount = 0;
  }

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*kstrst));

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)));
  long double testm = fabs(mKStarSigned) * mRSt * (1.0 - kstrst/(mRSt*fabs(mKStarSigned)));

  if ((testp> 15.0) && (testm> 15.0))
    {
      double fasymplus  = (1.0 - 1.0/((mRSt+kstrst)*pionac*mKStarSigned*mKStarSigned));
      double fasymminus = (1.0 - 1.0/((mRSt-kstrst)*pionac*mKStarSigned*mKStarSigned));
      return 0.5 * ((fasymplus + fasymminus)*cos(2*kstrst) + (2.0*sqrt(fasymplus*fasymminus)));
    }
  //    return (1.0 - 1.0/(mRSt*pionac*mKStarSigned*mKStarSigned))*(1.0+wavesign*cos(2*kstrst));
  
  dcomplex ffplus, ffminus;
  // Check for the classical limit in both functions separately
  if (((testp< 15.0) && (testm< 15.0))) // ||
    //       ((testp< 15.0) && (testm> 15.0) && (fabs(testp-testm < 1.0))) ||
    //       ((testp> 15.0) && (testm< 15.0) && (fabs(testp-testm < 1.0))))
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus);
      ccase = 1;
    }
  else if (testp< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus);
      //      GetFFsingle(&ffminus,-1);
      asym = (1.0 - 1.0/(mRSt*(1.0 - kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
      ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      ccase = 2;
    }
  else 
    {
      double asym;
      //      GetFFsingle(&ffplus);
      GetFFsingle(&ffminus, -1);
      asym = (1.0 - 1.0/(mRSt*(1.0 + kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      ccase = 3;
    }

  dcomplex expikr;
  expikr.re = cos(kstrst);
  expikr.im = sin(kstrst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));


  if (!finite(ffminus.re))
    cout << "FFMinus Re not a number !" << " " << ccase<< endl;
  
  if (!finite(ffminus.im))
    cout << "FFMinus Im not a number !" << " " << ccase<< endl;

  if ((ffplus.re > 2.0) || (ffplus.re < -2.0))
    cout << "FFplus Re wild !" << ffplus.re << endl;

  if ((ffplus.im > 2.0) || (ffplus.im < -2.0))
    cout << "FFplus Im wild !" << ffplus.im << endl;

  if ((ffminus.re > 2.0) || (ffminus.re < -2.0))
    cout << "FFminus Re wild !" << ffminus.re << " " << ccase << endl;
  
  if ((ffminus.im > 2.0) || (ffminus.im < -2.0))
    cout << "FFminus Im wild !" << ffminus.im << " " << ccase <<endl;

  coulqscpart = 0.5 * Gamow(fabs(mKStarSigned)) * (modl2(ffplus) + modl2(ffminus));

  return (0.5 * Gamow(fabs(mKStarSigned)) * 
	  (modl2(ffplus) + wavesign*sterm.re + wavesign*tterm.re + modl2(ffminus)));
}

// Calculate the wave-function modulus sqaured
// for non-identical particles
// that is the Coulomb interaction
double GetCoulomb()
{
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0);
  if (fabs(mRSt) == 0.0) return (Gamow(fabs(mKStarSigned)));

  // Classical limit - in the case of large k*r* product we go to 
  // classical coulomb interaction
  if (fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)))> 15.0)
    return (1.0 - 1.0/(mRSt*pionac*mKStarSigned*mKStarSigned));
  
  // Calculate the F function
  dcomplex ffplus;
  GetFFsingle(&ffplus);

  if (!finite(ffplus.re)) {
    cout << "FFPlus Re not a number !" << " " << endl;
    cout << mRSt << " " << mKStarSigned << " " << mROS << " " << mRSS << " " << mRLS << endl;
  }

  if (!finite(ffplus.im))
    cout << "FFPlus Im not a number !" << " " << endl;


  return (Gamow(fabs(mKStarSigned)) * modl2(ffplus));
}

// Calculate the wave-function modulus sqaured
// for non-identical particles
// for strong interaction in s-wave approximation
// without the Coulomb factor
double GetStrong()
{
  double fr, fi;
  double dr;
  double ar, ai;
  double ir, ii;
  
  double sr, si;
  double tr, ti;

  double skr, ckr;
  double srk, crk;

  double d0_sr = d0s.re;
  double d0_si = d0s.im;
  double d0_tr = d0t.re;
  double d0_ti = d0t.im;

  double f0_sr = f0s.re;
  double f0_si = f0s.im;
  double f0_tr = f0t.re;
  double f0_ti = f0t.im;

  double kstar = fabs(mKStarSigned);
  double rstar = mRSt;
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
/*
  if (pairtype == 8) {
    f0_sr = 2.88 /0.197327;
    f0_si = 0.0  /0.197327;
    d0_sr = 2.92 /0.197327;
    d0_si = 0.0  /0.197327;
    f0_tr = 1.66 /0.197327;
    f0_ti = 0.0  /0.197327;
    d0_tr = 3.78 /0.197327;
    d0_ti = 0.0  /0.197327;
  }
  else if (pairtype == 9) {
    f0_sr = -1.20 /0.197327;
    f0_si = 0.37  /0.197327;
    d0_sr = 0.0   /0.197327;
    d0_si = 0.0   /0.197327;
    f0_tr = -1.20 /0.197327;
    f0_ti = 0.37  /0.197327;
    d0_tr = 0.0   /0.197327;
    d0_ti = 0.0   /0.197327;
  }
*/
  if (rstar<1.0/0.197327) return 1.0;

  fr = f0_sr/(f0_sr*f0_sr + f0_si*f0_si);
  fi = -f0_si/(f0_sr*f0_sr + f0_si*f0_si);
  
  dr = 0.5*d0_sr*kstar*kstar;
  
  ir = fr+dr;
  ii = fi-kstar;
  
  ar = ir/(ir*ir+ii*ii);
  ai = -ii/(ir*ir+ii*ii);

  srk = sin(kstar*rstar);
  crk = cos(kstar*rstar);

  sr = (ar * crk - ai * srk) / rstar;
  si = (ai * crk + ar * srk) / rstar;
  
  fr = f0_tr/(f0_tr*f0_tr + f0_ti*f0_ti);
  fi = -f0_ti/(f0_tr*f0_tr + f0_ti*f0_ti);
  
  dr = 0.5*d0_tr*kstar*kstar;
  
  ir = fr+dr;
  ii = fi-kstar;
  
  ar = ir/(ir*ir+ii*ii);
  ai = -ii/(ir*ir+ii*ii);

  tr = (ar * crk - ai*srk) / rstar;
  ti = (ai * crk + ar*srk) / rstar;

  skr = sin(kstrst);
  ckr = cos(kstrst);

  return (0.25 * ((ckr+sr)*(ckr+sr) + (skr+si)*(skr+si)) + 
	  0.75 * ((ckr+tr)*(ckr+tr) + (skr+ti)*(skr+ti)));
}

// Calculates the h function for the strong interaction
long double geth(long double eta)
{
  long double etasum = log(1.0/eta) - euler;
  long double series = 0.0;
  long double x2inv = (eta*eta);
  long double element;
  long double save;
  for (int iter=1; iter<1000000; iter++) {
    element = ((1.0*iter)*(1.0*iter) + x2inv)*(1.0*iter);
    //    cout << "Element " << iter << " is " << element << endl;
    element = 1.0/element;
    if (iter == 1) save = 1.0e-10 * element;
    //    cout << "Element " << iter << " is " << element << endl;

    series += element;
    if (element < save) break;
  }
  series *= x2inv;
  etasum += series;

  return etasum;
}

long double chiim(long double eta)
{
  return Gamow(1.0/(eta*pionac))/(2.0*eta);
}

void bfunpfun(long double eta, long double rho, long double &bret, long double &pret)
{
  long double b0 = 1;
  long double b1 = eta*rho;
  long double bsum = b0+b1;
  long double bnpu, bn, bnmu;
  long double p0 = 1.0;
  long double p1 = 0.0;
  long double psum = p0+p1;
  long double pnpu, pn, pnmu;

  if (rho > FOURPI) {
    bret = sin(rho)/rho;
    pret = cos(rho);
    return;
  }
  

  bn = b1; bnmu = b0;
  pn = p1; pnmu = p0;
  for (int iter=1; iter<100000; iter++) {
    bnpu = 2 * eta*rho *bn - rho*rho*bnmu;
    bnpu /= (1.0*iter+1.0)*(1.0*iter+2.0);
    bsum += bnpu;
    //    cout << "B E " << iter << " " << bnpu << endl;

    pnpu = 2 * eta*rho *pn - rho*rho*pnmu - (2.0*iter+1.0)*2.0*eta*rho*bn;
    pnpu /= (1.0*iter)*(1.0*iter+1.0);
    psum += pnpu;
    //    cout << "P E " << iter << " " << pnpu << endl;

    bnmu = bn;
    bn   = bnpu;

    pnmu = pn;
    pn   = pnpu;
    if ((fabs(pnmu) + fabs (bnmu) + fabs(bnpu) + fabs(pnpu)) < 1.0e-20) 
      {
	//	cout << "iter " << iter << endl;
	break;
      }
    
  }

  

  bret = bsum;
  pret = psum;
}

// Calculate the Gtilde function
dcomplex GetG(long double eta, long double rho, long double hfun)
{
  dcomplex gtemp;
  long double bres, pres;
  long double bmult;

  bfunpfun(eta, rho, bres, pres);
  bmult = 2.0*eta*rho*bres;

  gtemp.re = pres + bmult*(log(fabs(2.0*eta*rho))+2.0*euler-1.0+hfun);
  gtemp.im = bmult * chiim(eta);

  return gtemp;
}

void Getfc(long double kstar, long double eta, long double hfun, dcomplex &fcs, dcomplex &fct)
{
  dcomplex ci;
  dcomplex cia;

  dcomplex fis;
  dcomplex dis;
  dcomplex fcinvs;

  dcomplex fit;
  dcomplex dit;
  dcomplex fcinvt;

  ci.re = hfun;
  ci.im = chiim(eta);
  cia = mult(ci, 2.0/pionac); 

  fis = invr(f0s);
  dis = mult(d0s, 0.5 * kstar * kstar);
  fcinvs.re = fis.re + dis.re - cia.re;
  fcinvs.im = fis.im + dis.im - cia.im;
  fcs = invr(fcinvs);
  
  fit = invr(f0t);
  dit = mult(d0t, 0.5 * kstar * kstar);
  fcinvt.re = fit.re + dit.re - cia.re;
  fcinvt.im = fit.im + dit.im - cia.im;
  fct = invr(fcinvt);
  
}


// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized) or fermions (antisymmetrized)
// with Coulomb interaction and Strong interaction included
double GetQuantumCoulombStrong()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  if (mRSt < 1.0/0.197327) 
    return GetQuantumCoulomb();
  
  double tKstRst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double kstar = fabs(mKStarSigned);
  long double rho = mRSt * kstar;
  
  int ccase = 0;
  static int pcount = 0;
  int wavesign = 1;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  //  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*tKstRst));

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = rho * (1.0 + tKstRst/(rho));
  long double testm = rho * (1.0 - tKstRst/(rho));

  dcomplex ffplus, ffminus;
  if ((testp> 15.0) && (testm> 15.0))
    {
      double asym;
      asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/rho)*pionac*kstar*kstar))/Gamow(kstar);
      //      cout << "as1 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
      ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);

      asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/rho)*pionac*kstar*kstar))/Gamow(kstar);
      //      cout << "as2 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      
    }
  
  // Check for the classical limit in both functions separately
  else if (((testp< 15.0) && (testm< 15.0))) // ||
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus);
      ccase = 1;
    }
  else if (testp< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus,1);
      GetFFsingle(&ffminus, -1);
      if ((fabs(ffminus.re) > 2.0) || fabs(ffminus.im) > 2.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/(rho)*pionac*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffminus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffminus.re = 1.0 + (asym -1.0) /2.0;
	ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      }
      ccase = 2;
    }
  else 
    {
      double asym;
      GetFFsingle(&ffminus, -1);
      GetFFsingle(&ffplus,   1);
      if ((fabs(ffplus.re) > 2.0) || fabs(ffplus.im) > 2.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/(rho)*pionac*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffplus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffplus.re = 1.0 + (asym -1.0) /2.0;
	ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      }
      ccase = 3;
    }

  long double eta  = 1.0/(kstar*pionac);
  long double hfun = geth(eta);
  dcomplex gtilde  = GetG(eta, rho, hfun);
  dcomplex gtilor  = mult(gtilde, 1.0/mRSt);

  dcomplex fcouls, fcoult;
  Getfc(kstar, eta, hfun, fcouls, fcoult);
  
  dcomplex fgs       = mult(gtilor, fcouls);
  long double fgmods = modl2(fgs);

  dcomplex fgt       = mult(gtilor, fcoult);
  long double fgmodt = modl2(fgt);

  dcomplex expikr;
  expikr.re = cos(tKstRst);
  expikr.im = -sin(tKstRst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));

  dcomplex epfpc = mult(expikr, ffplus);
  dcomplex emfmc = mult(expikrc, ffminus);

  long double fcgefhs = (fgs.re*emfmc.re + fgs.im*emfmc.im);
  long double fcgefgs = (fgs.re*epfpc.re + fgs.im*epfpc.im);

  long double fcgefht = (fgt.re*emfmc.re + fgt.im*emfmc.im);
  long double fcgefgt = (fgt.re*epfpc.re + fgt.im*epfpc.im);

  long double smult = 1+wavesign;

  if (!finite(ffminus.re))
    cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase<< endl;
  
  if (!finite(ffminus.im))
    cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase<< endl;

  if ((ffplus.re > 2.0) || (ffplus.re < -2.0))
    cout << "FFplus Re wild !" << ffplus.re << " case " << ccase << " " << testp << " " << testm << endl;

  if ((ffplus.im > 2.0) || (ffplus.im < -2.0))
    cout << "FFplus Im wild !" << ffplus.im << " case " << ccase << " " << testp << " " << testm << endl;

  if ((ffminus.re > 2.0) || (ffminus.re < -2.0))
    cout << "FFminus Re wild !" << ffminus.re << " case " << ccase << endl;
  
  if ((ffminus.im > 2.0) || (ffminus.im < -2.0))
    cout << "FFminus Im wild !" << ffminus.im << " case " << ccase <<endl;

  //  coulqscpart = 0.5 * Gamow(kstar) * (modl2(ffplus) + modl2(ffminus));

  //   return (0.5 * Gamow(kstar) * 
  // 	  (modl2(ffplus) + wavesign*sterm.re + wavesign*tterm.re + modl2(ffminus)));
  if (twospin == 1) {
    wavesign = 1;
    smult = 2;
    long double singlet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefhs + fcgefgs)));
    wavesign = -1;
    smult = 0;
    long double triplet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmodt * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefht + fcgefgt)));
    //    scmp = singlet;
    //    tcmp = triplet;

    //    ccmp = 0.5 * Gamow(kstar) * (modl2(ffplus) + modl2(ffminus));
    //    gcmp = fgmod;

    return (0.25 * singlet + 0.75 * triplet);
    //    return triplet;
  }
  else
    return (0.5 * Gamow(kstar) *
	    (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
	     wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefhs + fcgefgs)));
}

// Calculate the wave-function modulus sqaured
// with Coulomb interaction and Strong interaction included
double GetCoulombStrong()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  if (mRSt < 1.0/0.197327) 
    return GetCoulomb();
  
  double tKstRst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double kstar = fabs(mKStarSigned);
  long double rho = mRSt * kstar;
  
  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = rho * (1.0 + tKstRst/(rho));

  dcomplex ffplus;
  if (testp> 15.0)
    {
      double asym;
      asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/rho)*pionac*kstar*kstar))/Gamow(kstar);
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
    }
  else 
    GetFFsingle(&ffplus,1);

  
  long double eta  = 1.0/(kstar*pionac);
  long double hfun = geth(eta);
  dcomplex gtilde  = GetG(eta, rho, hfun);
  dcomplex gtilor  = mult(gtilde, 1.0/mRSt);

  dcomplex fcouls, fcoult;
  Getfc(kstar, eta, hfun, fcouls, fcoult);
  
  dcomplex fgs       = mult(gtilor, fcouls);
  //  long double fgmods = modl2(fgs);

  dcomplex fgt       = mult(gtilor, fcoult);
  //  long double fgmodt = modl2(fgt);

  dcomplex expikr;
  expikr.re = cos(tKstRst);
  expikr.im = -sin(tKstRst);
  //  dcomplex expikrc = conj(expikr);
  //  dcomplex ffplusc = conj(ffplus);

  //  dcomplex expikr2 = mult(expikr, expikr);
  //  dcomplex expikrc2 = conj(expikr2);

  dcomplex epfpc = mult(expikr, ffplus);
  //  dcomplex emfmc = mult(expikrc, ffminus);

//   long double fcgefhs = (fgs.re*emfmc.re + fgs.im*emfmc.im);
//  long double fcgefgs = (fgs.re*epfpc.re + fgs.im*epfpc.im);
//  long double fcgefgt = (fgt.re*epfpc.re + fgt.im*epfpc.im);
  //  long double fcgefht = (fgt.re*emfmc.re + fgt.im*emfmc.im);

  dcomplex fvfs;
  fvfs.re = epfpc.re + fgs.re;
  fvfs.im = epfpc.im + fgs.im;

  dcomplex fvft;
  fvft.re = epfpc.re + fgt.re;
  fvft.im = epfpc.im + fgt.im;

//   long double smult = 1+wavesign;

//  return 0.5 * Gamow(kstar);
//   return (0.5 * Gamow(kstar) *
// 	  (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
// 	   wavesign*sterm.re + wavesign*tterm.re +
// 	   smult * 2 * (fcgefhs + fcgefgs)));

  if (twospin == 1) {
  //   wavesign = 1;
//     smult = 2;
//     long double singlet = (0.5 * Gamow(kstar) *
// 			   (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
// 			    wavesign*sterm.re + wavesign*tterm.re +
// 			    smult * 2 * (fcgefhs + fcgefgs)));
//     wavesign = -1;
//     smult = 0;
//     long double triplet = (0.5 * Gamow(kstar) *
// 			   (2.0 * fgmodt * smult + modl2(ffplus) + modl2(ffminus) + 
// 			    wavesign*sterm.re + wavesign*tterm.re +
// 			    smult * 2 * (fcgefht + fcgefgt)));
    //    scmp = singlet;
    //    tcmp = triplet;

    //    ccmp = 0.5 * Gamow(kstar) * (modl2(ffplus) + modl2(ffminus));
    //    gcmp = fgmod;

    long double singlet = Gamow(kstar) * modl2(fvfs);
    long double triplet = Gamow(kstar) * modl2(fvft);

    //    cout << " s t " << singlet << "  " << triplet << "   -   "  << (0.25*singlet+0.75*triplet) << endl;

    return (0.25 * singlet + 0.75 * triplet);
    //    return triplet;
  }
  else
    return Gamow(kstar) * modl2(fvfs);
//     return (0.5 * Gamow(kstar) *
// 	    (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
// 	     wavesign*sterm.re + wavesign*tterm.re +
// 			    smult * 2 * (fcgefhs + fcgefgs)));

}


double 
funeh(double xarg, double rad, double alfa)
{
  return exp(-sqrt(xarg*xarg/(rad*rad) + alfa*alfa));
}

double 
funex(double xarg, double rad)
{
  return exp(-xarg/rad);
  
}

//void PairKinematics(PARTICLE part1, PARTICLE part2)
void PairKinematics(ParticleCoor part1, ParticleCoor part2)
{
  if (mMomRes1 && mMomRes2) {
    // Calculate pair variables
    // Taking into accout STAR's momentum resolution
    // First calculate the particles momenta after smearing

    double tPhi = atan2(part1.py,part1.px);
    double tP = part1.px*part1.px+part1.py*part1.py+part1.pz*part1.pz;
    tP = sqrt(tP);

    double per = mMomRes1->PResA+mMomRes1->PResB*pow(tP,mMomRes1->PResAlpha)+mMomRes1->PResC*tP;
    double thetaan = TMath::ATan2(hypot(part1.px,part1.py),part1.pz);
    double phier = mMomRes1->PhiA+mMomRes1->PhiB*pow(tP,mMomRes1->PhiAlpha);
    double thetaer = mMomRes1->ThetaA+mMomRes1->ThetaB*pow(tP,mMomRes1->ThetaAlpha);
    double pshift = (mMomRes1->PMeanA+mMomRes1->PMeanB*pow(tP,mMomRes1->PMeanAlpha));
    double tmass = TMath::Sqrt(part1.e*part1.e - tP*tP);

    
    // Apply the momentum resolution in reverse
    // Here we have REAL momenta and want to "truify" them
    double ptrue = tP + pshift; // * mMomResScale;
    
    double ptruex = ptrue * part1.px / tP;
    double ptruey = ptrue * part1.py / tP;
    double ptruez = ptrue * part1.pz / tP;

    double Deltapx = TMath::Abs(ptruex) * per + TMath::Abs(ptruey) * phier + TMath::Abs(ptruex * (1/TMath::Tan(thetaan))) * thetaer;
    double Deltapy = TMath::Abs(ptruey) * per + TMath::Abs(ptruex) * phier + TMath::Abs(ptruey * (1/TMath::Tan(thetaan))) * thetaer;
    double Deltapz = TMath::Abs(ptruez) * per + TMath::Abs(ptruez * TMath::Tan(thetaan)) * thetaer;

    double p1x = ((ptruex + mRand->Gaus(0,fabs(Deltapx) * mMomResScale)));
    double p1y = ((ptruey + mRand->Gaus(0,fabs(Deltapy) * mMomResScale)));
    double p1z = ((ptruez + mRand->Gaus(0,fabs(Deltapz) * mMomResScale)));
    double p1e = TMath::Sqrt(p1x*p1x +
 			     p1y*p1y +
 			     p1z*p1z +
 			     tmass*tmass);
    
    tPhi = atan2(part2.py,part2.px);
    tP = part2.px*part2.px+part2.py*part2.py+part2.pz*part2.pz;
    tP = sqrt(tP);
    
    per = mMomRes2->PResA+mMomRes2->PResB*pow(tP,mMomRes2->PResAlpha)+mMomRes2->PResC*tP;
    thetaan = TMath::ATan2(hypot(part2.px,part2.py),part2.pz);
    phier = mMomRes2->PhiA+mMomRes2->PhiB*pow(tP,mMomRes2->PhiAlpha);
    thetaer = mMomRes2->ThetaA+mMomRes2->ThetaB*pow(tP,mMomRes2->ThetaAlpha);
    pshift = (mMomRes2->PMeanA+mMomRes2->PMeanB*pow(tP,mMomRes2->PMeanAlpha));
    tmass = TMath::Sqrt(part2.e*part2.e - tP*tP);
    
    ptrue = tP + pshift;// * mMomResScale;
    
    ptruex = ptrue * part2.px / tP;
    ptruey = ptrue * part2.py / tP;
    ptruez = ptrue * part2.pz / tP;
    
    Deltapx = TMath::Abs(ptruex) * per + TMath::Abs(ptruey) * phier + TMath::Abs(ptruex * (1/TMath::Tan(thetaan))) * thetaer;
    Deltapy = TMath::Abs(ptruey) * per + TMath::Abs(ptruex) * phier + TMath::Abs(ptruey * (1/TMath::Tan(thetaan))) * thetaer;
    Deltapz = TMath::Abs(ptruez) * per + TMath::Abs(ptruez * TMath::Tan(thetaan)) * thetaer;
    
    double p2x = ((ptruex + mRand->Gaus(0,fabs(Deltapx) * mMomResScale)));
    double p2y = ((ptruey + mRand->Gaus(0,fabs(Deltapy) * mMomResScale)));
    double p2z = ((ptruez + mRand->Gaus(0,fabs(Deltapz) * mMomResScale)));
    double p2e = TMath::Sqrt(p2x*p2x +
			     p2y*p2y +
			     p2z*p2z +
			     tmass*tmass);
    
    double tPxm = p1x+p2x;
    double tPym = p1y+p2y;
    double tPzm = p1z+p2z;
    double tEm  = p1e+p2e;
    double tPtm = tPxm*tPxm + tPym*tPym;
    double tMtm = tEm*tEm - tPzm*tPzm;//mCVK;
    double tMm  = sqrt(tMtm - tPtm);
    tMtm = sqrt(tMtm);
    tPtm = sqrt(tPtm);

    mKTm = tPtm/2.0;
    pairphim = TMath::ATan2(tPym, tPxm);

    if ((mKTm<ktmin) || (mKTm>ktmax))
      return;
    
    mBetatm = tPtm/tMtm;
    
    if ((mBetatm<btmin) || (mBetatm>btmax)) {
      mKTm = -1.0;
      return;
    }
    
    // Boost to LCMS
    double tBetam = tPzm/tEm;
    double tGammam = tEm/tMtm;	    
    mKStarLongm = tGammam * (p1z - tBetam * p1e);
    double tE1Lm = tGammam * (p1e  - tBetam * p1z);
    
    particle1prfm.pz = particle1lcmsm.pz = tGammam * (p1z - tBetam * p1e);
    particle1lcmsm.e  = tGammam * (p1e  - tBetam * p1z);
  
    particle2prfm.pz = particle2lcmsm.pz = tGammam * (p2z - tBetam * p2e);
    particle2lcmsm.e  = tGammam * (p2e  - tBetam * p2z);
  
    // Rotate in transverse plane
    mKStarOutm  = ( p1x*tPxm + p1y*tPym)/tPtm;
    mKStarSidem = (-p1x*tPym + p1y*tPxm)/tPtm;
      
    particle1lcmsm.px = mKStarOutm;
    particle1lcmsm.py = mKStarSidem;
  
    particle2lcmsm.px = (p2x*tPxm + p2y*tPym)/tPtm;
    particle2lcmsm.py = (p2y*tPxm - p2x*tPym)/tPtm;
  
    mKOm = particle1lcmsm.px - particle2lcmsm.px;
    mKSm = particle1lcmsm.py - particle2lcmsm.py;
    mKLm = particle1lcmsm.pz - particle2lcmsm.pz;
    mDEm = particle1lcmsm.e  - particle2lcmsm.e;
  
    // Boost to pair cms
    mKStarOutm = tMtm/tMm * (mKStarOutm - tPtm/tMtm * tE1Lm);
  
    Double_t tBetatm = tPtm/tMtm;
    Double_t tGammatm = 1.0/TMath::Sqrt(1.0-tBetatm*tBetatm);
  
    particle1prfm.x = tGammatm*(particle1lcmsm.x - tBetatm*particle1lcmsm.t);
    particle1prfm.t = tGammatm*(particle1lcmsm.t - tBetatm*particle1lcmsm.x);

    particle2prfm.x = tGammatm*(particle2lcmsm.x - tBetatm*particle2lcmsm.t);
    particle2prfm.t = tGammatm*(particle2lcmsm.t - tBetatm*particle2lcmsm.x);
  
    mKOm = particle1lcmsm.px - particle2lcmsm.px;

    mKStarSignedm = mKStarOutm>0.? 1. : -1.;
    mKStarSignedm *= sqrt(mKStarSidem*mKStarSidem + mKStarOutm*mKStarOutm + mKStarLongm*mKStarLongm);
    
  }
  // Calculate pair variables
  double tPx = part1.px+part2.px;
  double tPy = part1.py+part2.py;
  double tPz = part1.pz+part2.pz;
  double tE  = part1.e +part2.e;
  double tPt = tPx*tPx + tPy*tPy;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  mKT = tPt/2.0;
  pairphi = TMath::ATan2(tPy, tPx);

  if (tE<0.0) { 
    cout << "Negative total energy !!! " << endl; 
    cout << part1.px << " " << part1.py << " " << part1.pz << " " << part1.e << endl;
    cout << part2.px << " " << part2.py << " " << part2.pz << " " << part2.e << endl;
    exit(0); 
  }

  mBetat = tPt/tMt;
    
  //   if ((mBetat > btmin) && (mBetat < btmax))
  //     {
  // Boost to LCMS
  double tBeta = tPz/tE;
  double tGamma = tE/tMt;	    
  mKStarLong = tGamma * (part1.pz - tBeta * part1.e);
  double tE1L = tGamma * (part1.e  - tBeta * part1.pz);
  
  // Transform to LCMS
  particle1lcms.z = tGamma * (part1.z - tBeta * part1.t);
  particle1lcms.t = tGamma * (part1.t - tBeta * part1.z);
  
  particle2lcms.z = tGamma * (part2.z - tBeta * part2.t);
  particle2lcms.t = tGamma * (part2.t - tBeta * part2.z);
  
  particle1prf.pz = particle1lcms.pz = tGamma * (part1.pz - tBeta * part1.e);
  particle1lcms.e  = tGamma * (part1.e  - tBeta * part1.pz);
  
  particle2prf.pz = particle2lcms.pz = tGamma * (part2.pz - tBeta * part2.e);
  particle2lcms.e  = tGamma * (part2.e  - tBeta * part2.pz);
  
  // Rotate in transverse plane
  mKStarOut  = ( part1.px*tPx + part1.py*tPy)/tPt;
  mKStarSide = (-part1.px*tPy + part1.py*tPx)/tPt;
      
  particle1lcms.px = mKStarOut;
  particle1lcms.py = mKStarSide;
  
  particle2lcms.px = (part2.px*tPx + part2.py*tPy)/tPt;
  particle2lcms.py = (part2.py*tPx - part2.px*tPy)/tPt;;
  
  mKO = particle1lcms.px - particle2lcms.px;
  mKS = particle1lcms.py - particle2lcms.py;
  mKL = particle1lcms.pz - particle2lcms.pz;
  mDE = particle1lcms.e  - particle2lcms.e;
  
  // save the rotated coordinates in LCMS variables
  particle1lcms.x = ( part1.x*tPx + part1.y*tPy)/tPt;
  particle1lcms.y = (-part1.x*tPy + part1.y*tPx)/tPt;

  particle2lcms.x = ( part2.x*tPx + part2.y*tPy)/tPt;
  particle2lcms.y = (-part2.x*tPy + part2.y*tPx)/tPt;
  
  // Boost to pair cms
  mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);
  
  Double_t tBetat = tPt/tMt;
  Double_t tGammat = 1.0/TMath::Sqrt(1.0-tBetat*tBetat);
  
  particle1prf.x = tGammat*(particle1lcms.x - tBetat*particle1lcms.t);
  particle1prf.t = tGammat*(particle1lcms.t - tBetat*particle1lcms.x);
  
  particle2prf.x = tGammat*(particle2lcms.x - tBetat*particle2lcms.t);
  particle2prf.t = tGammat*(particle2lcms.t - tBetat*particle2lcms.x);
  
  mRO = (particle1lcms.x - particle2lcms.x)/0.197327;
  mRS = (particle1lcms.y - particle2lcms.y)/0.197327;
  mRL = (particle1lcms.z - particle2lcms.z)/0.197327;
  mDT = (particle1lcms.t - particle2lcms.t)/0.197327;

  mKStarSigned = mKStarOut>0.? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide*mKStarSide + mKStarOut*mKStarOut + mKStarLong*mKStarLong);

  double tDX = part1.x-part2.x;
  double tDY = part1.y-part2.y;
  mRLong = part1.z-part2.z;
  mDTime = part1.t-part2.t;

  mRTrans = tDX>0.? ::sqrt(tDX*tDX+tDY*tDY) : -1.*::sqrt(tDX*tDX+tDY*tDY);
  mROut = (tDX*tPx + tDY*tPy)/tPt;
  mRSide = (-tDX*tPy + tDY*tPx)/tPt;

  mRSidePairCMS = mRSide;
  mRSS = mRSidePairCMS/0.197327;

  mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
  mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);

  mRLS = mRLongPairCMS/0.197327;
  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
  mROS = mROutPairCMS/0.197327;
  mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);

  mRStar = ::sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS +
		  mRLongPairCMS*mRLongPairCMS);
  mRSt = mRStar/0.197327;

  ev2weight = part1.eventweight * part2.eventweight * 10e10;

  Double_t tPhi1 = atan2(part1.py,part1.px);
  Double_t tPhi2 = atan2(part2.py,part2.px);

  Double_t tEta1 = 0.5*log((part1.e+part1.pz)/(part1.e-part1.pz));
  Double_t tEta2 = 0.5*log((part2.e+part2.pz)/(part2.e-part2.pz));

  deta = tEta2 - tEta1;
  dphi = tPhi2 - tPhi1;

  if (dphi <  phiLow) dphi += TMath::Pi()*2.0;
  if (dphi > phiUp) dphi -= TMath::Pi()*2.0;

  if ((!mMomRes1) || (!mMomRes2)) {
    mKTm = mKT;
    pairphim = pairphi;

    if ((mKTm<ktmin) || (mKTm>ktmax))
      return;
    
    mBetatm = mBetat;
    
    if ((mBetatm<btmin) || (mBetatm>btmax)) {
      mKTm = -1.0;
      return;
    }
    
    mKStarLongm = mKStarLong;
    
    particle1prfm.pz = particle1prf.pz;
    particle1lcmsm.e  = particle1lcms.e;
  
    particle2prfm.pz = particle2prf.pz;
    particle2lcmsm.e  = particle2lcms.e;
  
    mKStarOutm  = mKStarOut;
    mKStarSidem = mKStarSide;
      
    particle1lcmsm.px = mKStarOutm;
    particle1lcmsm.py = mKStarSidem;
  
    particle2lcmsm.px = particle2lcms.px;
    particle2lcmsm.py = particle2lcms.py;
  
    mKOm = mKO;
    mKSm = mKS;
    mKLm = mKL;
    mDEm = mDE;
  
    mKStarOutm = mKStarOut;
  
    particle1prfm.x = particle1prf.x;
    particle1prfm.t = particle1prf.t;
    			      				  
    particle2prfm.x = particle2prf.x;
    particle2prfm.t = particle2prf.t;
  
    //    mKOm = particle1lcmsm.px - particle2lcmsm.px;

    mKStarSignedm = mKStarOutm>0.? 1. : -1.;
    mKStarSignedm *= sqrt(mKStarSidem*mKStarSidem + mKStarOutm*mKStarOutm + mKStarLongm*mKStarLongm);
    
  }

  //TPISourceMonitorFill();
}

double GetQuantum()
{
  double quantumweight;

  if(partpid2 != 0 && partpid!=partpid2) return 1;

  if (twospin == 0) {
    quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
  }
  else if (twospin == 1) {
    if (pcount ==3) {
      quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
    }
    else {
      quantumweight = 1.0-TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
    }
    pcount++;
    if (pcount == 4) pcount=0;
  }
//   if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
//     quantumweight = 1.0;

  return quantumweight;
}

double GetFull()
{
  double coulombweight;
//   if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0) {
//     return 1.0;
//   }
  if (partpid2 != 0) {
    if (fabs(pionac) > 0.1) {
//       if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
// 	coulombweight = 1.0;
//       else
      
      if (fabs(d0s.re) > 0.00001)
	coulombweight = GetCoulombStrong();
      else
	coulombweight = GetCoulomb();
    }
    else {
      coulombweight = GetStrong();
    }
  }
  else {
    if (pairtype != 12)
      coulombweight = GetQuantumCoulomb();
    else 
      coulombweight = GetQuantumCoulombStrong();
  }

//   if (coulombweight > 10.0) {
//     cout << "Weigth for " << fabs(mKStarSigned) << " " << fiter << " " << siter << " " << " is " << coulombweight  << endl;
//   }
			
  return coulombweight;
}

double GetFullNoQS()
{
  double coulombweight;
//   if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0) {
//     return 1.0;
//   }
  if (partpid2 != 0) {
    if (fabs(pionac) > 0.1) {
//       if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
// 	coulombweight = 1.0;
//       else
      
      if (fabs(d0s.re) > 0.00001)
	coulombweight = GetCoulombStrong();
      else
	coulombweight = GetCoulomb();
    }
    else {
      coulombweight = GetStrong();
    }
  }
  else {
    if (pairtype != 12)
      coulombweight = GetCoulomb();
    else 
      coulombweight = GetCoulombStrong();
  }

//   if (coulombweight > 10.0) {
//     cout << "Weigth for " << fabs(mKStarSigned) << " " << fiter << " " << siter << " " << " is " << coulombweight  << endl;
//   }
			
  return coulombweight;
}
