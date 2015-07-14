#include "SourceMonitor.h"

//extern PARTICLE particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf, buf;
extern ParticleCoor particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf, buf;
extern double mRL, mRS, mRO, mDT;
extern double mKO, mKS, mKL, mDE;
extern double mKOm, mKSm, mKLm, mDEm;
extern double mROS, mRSS, mRLS, mRSt;
extern double mKStarLong, mKStarOut, mKStarSide, mKStarSigned, mRLong, mDTime, mRTrans, mROut, mRSide, mRSidePairCMS, mRLongPairCMS, mDTimePairLCMS, mROutPairCMS, mDTimePairCMS, mRStar, mBetat;
extern double mKStarLongm, mKStarOutm, mKStarSidem, mKStarSignedm, mBetatm;
extern double pairphi;
extern double pairphim;
extern double pairweight;

TH1D *tROutLCMS;
TH1D *tRSideLCMS;
TH1D *tRLongLCMS;
TH1D *tRTimeLCMS;

TH1D *tROutPRF;
TH1D *tRSidePRF;
TH1D *tRLongPRF;
TH1D *tRTimePRF;

TH3D *tROutSideLongPRF;
TH3D *tROutSideLongLCMS;

TH2D *tROutLCMSQinv;
TH2D *tRSideLCMSQinv;
TH2D *tRLongLCMSQinv;
TH2D *tRTimeLCMSQinv;

TH2D *tROutPRFQinv;
TH2D *tRSidePRFQinv;
TH2D *tRLongPRFQinv;
TH2D *tRTimePRFQinv;

TH1D *tRInvPRF;

TH2D *simage1;
TH2D *simage2;

TH1D *time1;
TH1D *time2;

TH2D *tMomResX;
TH2D *tMomResY;
TH2D *tMomResZ;

Bool_t isoff;

void TPISourceMonitorInit()
{
  isoff = 0;
  
  tROutLCMS = new TH1D("ROutLCMS", "ROut in LCMS", 400, -20.0, 20.0);
  tRSideLCMS = new TH1D("RSideLCMS", "RSide in LCMS", 400, -20.0, 20.0);
  tRLongLCMS = new TH1D("RLongLCMS", "RLong in LCMS", 400, -20.0, 20.0);
  tRTimeLCMS = new TH1D("RTimeLCMS", "DTime in LCMS", 400, -20.0, 20.0);
  
  tROutSideLongLCMS = new TH3D("ROutSideLongLCMS", "Space distribution in LCMS", 200, -100.0, 100.0, 200, -100.0, 100.0, 200, -100.0, 100.0);

  tROutPRF = new TH1D("ROutPRF", "ROut in PRF", 400, -20.0, 20.0);
  tRSidePRF = new TH1D("RSidePRF", "RSide in PRF", 400, -20.0, 20.0);
  tRLongPRF = new TH1D("RLongPRF", "RLong in PRF", 400, -20.0, 20.0);
  tRTimePRF = new TH1D("RTimePRF", "DTime in PRF", 400, -20.0, 20.0);

  tROutSideLongPRF = new TH3D("ROutSideLongPRF", "Space distribution in PRF", 200, -100.0, 100.0, 200, -100.0, 100.0, 200, -100.0, 100.0);

  tROutLCMSQinv  = new TH2D("ROutLCMSKStar",  "ROut in LCMS vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRSideLCMSQinv = new TH2D("RSideLCMSKStar", "RSide in LCMS vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRLongLCMSQinv = new TH2D("RLongLCMSKStar", "RLong in LCMS vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRTimeLCMSQinv = new TH2D("RTimeLCMSKStar", "DTime in LCMS vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);

  tROutPRFQinv  = new TH2D("ROutPRFKStar",  "ROut in PRF vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRSidePRFQinv = new TH2D("RSidePRFKStar", "RSide in PRF vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRLongPRFQinv = new TH2D("RLongPRFKStar", "RLong in PRF vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);
  tRTimePRFQinv = new TH2D("RTimePRFKStar", "DTime in PRF vs qinv", 400, -20.0, 20.0, 60, 0.0, 0.3);

  tRInvPRF = new TH1D("RInvPRF", "RInvariant", 400, -20.0, 20.0);

  simage1 = new TH2D("simage1", "simage1", 100, -15.0, 15.0, 100, -15.0, 15.0);
  simage2 = new TH2D("simage2", "simage2", 100, -15.0, 15.0, 100, -15.0, 15.0);

  time1 = new TH1D("time1", "Emission time particle1;t [fm/c]; dN/dt", 100, 0.0, 50.0);
  time2 = new TH1D("time2", "Emission time particle2;t [fm/c]; dN/dt", 100, 0.0, 50.0);

  tMomResX = new TH2D("MomResX","Momentum resolution out", 50, 0.0, 0.2, 100, -0.1, 0.1);
  tMomResY = new TH2D("MomResY","Momentum resolution side", 50, 0.0, 0.2, 100, -0.1, 0.1);
  tMomResZ = new TH2D("MomResZ","Momentum resolution long", 50, 0.0, 0.2, 100, -0.1, 0.1);

  tROutLCMS->Sumw2();
  tRSideLCMS->Sumw2();
  tRLongLCMS->Sumw2();
  tRTimeLCMS->Sumw2();

  tROutSideLongLCMS->Sumw2();

  tROutPRF->Sumw2();
  tRSidePRF->Sumw2();
  tRLongPRF->Sumw2();
  tRTimePRF->Sumw2();

  tROutSideLongPRF->Sumw2();

  tROutLCMSQinv->Sumw2();
  tRSideLCMSQinv->Sumw2();
  tRLongLCMSQinv->Sumw2();
  tRTimeLCMSQinv->Sumw2();

  tROutPRFQinv->Sumw2();
  tRSidePRFQinv->Sumw2();
  tRLongPRFQinv->Sumw2();
  tRTimePRFQinv->Sumw2();

  tRInvPRF->Sumw2();
}

void TPISourceMonitorFill()
{
  if (!isoff) {
    tROutLCMS ->Fill(mROut);
    tRSideLCMS->Fill(mRSide);
    tRLongLCMS->Fill(mRLong);
    tRTimeLCMS->Fill(mDTimePairLCMS);
    tROutSideLongLCMS->Fill(mROut, mRSide, mRLong);
  
    tROutPRF ->Fill(mROutPairCMS);
    tRSidePRF->Fill(mRSidePairCMS);
    tRLongPRF->Fill(mRLongPairCMS);
    tRTimePRF->Fill(mDTimePairCMS);
    tROutSideLongPRF->Fill(mROutPairCMS, mRSidePairCMS, mRLongPairCMS);

    tROutLCMSQinv->Fill(mROut, fabs(mKStarSigned)*2.0);
    tRSideLCMSQinv->Fill(mRSide, fabs(mKStarSigned)*2.0);
    tRLongLCMSQinv->Fill(mRLong, fabs(mKStarSigned)*2.0);
    tRTimeLCMSQinv->Fill(mDTimePairLCMS, fabs(mKStarSigned)*2.0);

    tROutPRFQinv->Fill(mROutPairCMS, fabs(mKStarSigned)*2.0);
    tRSidePRFQinv->Fill(mRSidePairCMS, fabs(mKStarSigned)*2.0);
    tRLongPRFQinv->Fill(mRLongPairCMS, fabs(mKStarSigned)*2.0);
    tRTimePRFQinv->Fill(mDTimePairCMS, fabs(mKStarSigned)*2.0);

    tRInvPRF->Fill(mRStar, 1.0/(mRStar*mRStar));  

    if (fabs(mKStarSignedm) < 0.1) {
      double rad = hypot(particle1lcms.x, particle1lcms.y);
      double phi = TMath::ATan2(particle1lcms.y, particle1lcms.x);
      simage1->Fill(rad * cos(phi - pairphi), rad * sin(phi - pairphi), 1.0);
      time1->Fill(particle1lcms.t);
    
      rad = hypot(particle2lcms.x, particle2lcms.y);
      phi = TMath::ATan2(particle2lcms.y, particle2lcms.x);
      simage2->Fill(rad * cos(phi - pairphi), rad * sin(phi - pairphi), 1.0);
      time2->Fill(particle2lcms.t);
    }
    
    tMomResX->Fill(fabs(mKStarSigned), mKStarOut - mKStarOutm);
    tMomResY->Fill(fabs(mKStarSigned), mKStarSide - mKStarSidem);
    tMomResZ->Fill(fabs(mKStarSigned), mKStarLong - mKStarLongm);
  }
}

void TPISourceMonitorWrite()
{
  if (!isoff) {
    tROutLCMS->Write();
    tRSideLCMS->Write();
    tRLongLCMS->Write();
    tRTimeLCMS->Write();
    
    tROutSideLongLCMS->Write();
    
    tROutPRF->Write();
    tRSidePRF->Write();
    tRLongPRF->Write();
    tRTimePRF->Write();
    
    tROutSideLongPRF->Write();
    
    tROutLCMSQinv->Write();
    tRSideLCMSQinv->Write();
    tRLongLCMSQinv->Write();
    tRTimeLCMSQinv->Write();
    
    tROutPRFQinv->Write();
    tRSidePRFQinv->Write();
    tRLongPRFQinv->Write();
    tRTimePRFQinv->Write();
    
    tRInvPRF->Write();
    simage1->Write();
    simage2->Write();
    
    time1->Write();
    time2->Write();

    tMomResX->Write();
    tMomResY->Write();
    tMomResZ->Write();
  }
}

void TPISourceMonitorSetOff(Bool_t dooff)
{
  isoff = dooff;
}
