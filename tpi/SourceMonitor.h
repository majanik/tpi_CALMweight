#ifndef _TPI_SourceMonitor_
#define _TPI_SourceMonitor_

#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>

#include "TPIGlobal.h"
#include "PairWeight.h"

void TPISourceMonitorInit();
void TPISourceMonitorFill();
void TPISourceMonitorWrite();
void TPISourceMonitorSetOff(Bool_t dooff);

#endif
