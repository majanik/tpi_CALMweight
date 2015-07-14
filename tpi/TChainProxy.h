#ifndef _TCHAINPROXY_H_
#define _TCHAINPROXY_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include "ParticleCoor.h"

const unsigned int _maxParticleCount = 1e5;

class TChainProxy
{
private:
	enum ChainType
	{
		Therminator = 0,
		Siniukow = 1,
		EPOS = 2,
		NONE = 99
	};


	const char *_types[3]; // TTree names in files
	std::vector<Int_t> _particlesPerEvent;
	Int_t _iEvent;
	Int_t _totalEventCount;

	ChainType _type;
	TChain *_chain;

	ParticleCoor *_therminatorParticle;

	/* 
	 * Siniukow
	 */
	typedef struct
	{
		Int_t		npart;
		Int_t		id[_maxParticleCount];
		Int_t		mid[_maxParticleCount];
		Float_t		x[_maxParticleCount];
		Float_t		y[_maxParticleCount];
		Float_t		z[_maxParticleCount];
		Float_t		t[_maxParticleCount];
		Float_t		px[_maxParticleCount];
		Float_t		py[_maxParticleCount];
		Float_t		pz[_maxParticleCount];
		Float_t		E[_maxParticleCount];
	}SiniukowEvent;
	std::vector<ParticleCoor> _siniukowParticles;
	Bool_t _isSiniukowLoaded;
	void _loadSiniukow();

	/*
	 * EPOS
	 */
	typedef struct
	{
		Int_t		np;
		Float_t		bim;
		Float_t		zus[_maxParticleCount];
		Float_t		px[_maxParticleCount];
		Float_t		py[_maxParticleCount];
		Float_t		pz[_maxParticleCount];
		Float_t		e[_maxParticleCount];
		Float_t		x[_maxParticleCount];
		Float_t		y[_maxParticleCount];
		Float_t		z[_maxParticleCount];
		Float_t 	t[_maxParticleCount];
		Int_t		id[_maxParticleCount];
		Int_t		ist[_maxParticleCount];
		Int_t 		ity[_maxParticleCount];
		Int_t 		ior[_maxParticleCount];
		Int_t 		jor[_maxParticleCount];
	}EPOSEvent;
	const int _EPOSeventsPerFile;
	std::vector<ParticleCoor> _EPOSParticles;
	Bool_t _isEPOSLoaded;
	std::vector<unsigned int> _EPOSeventIndices;
	void _loadEPOS();
	void _loadEPOSfilesChain();
	Int_t _translateEPOSpidToTherminator(Int_t pid);

public:
	TChainProxy();
	~TChainProxy();
	Int_t Add(const char *fileName);
	
	// ALWAYS RETURNS 1 !
	Int_t GetEntry(Long64_t entryNumber);

	// ALWAYS RETURNS 1 !
	Int_t SetBranchAddress(const char *skip, ParticleCoor *particle);

	Long64_t GetEntries();

};

#endif