#ifndef _TPI_PairWeight_
#define _TPI_PairWeight_

#include <iostream>
#include "TPIGlobal.h"

using namespace std;

void InitializeGamow();
//void PairKinematics(PARTICLE part1, PARTICLE part2);
void PairKinematics(ParticleCoor part1, ParticleCoor part2);

double GetQuantum();
double GetFull();
double GetFullNoQS();

#endif
