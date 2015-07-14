#include "TChainProxy.h"

TChainProxy::TChainProxy() : _type(NONE), _isSiniukowLoaded(false), _iEvent(0), _EPOSeventsPerFile(12)
{
	_types[0]  = "particles";	// Therminator
	_types[1]  = "tf";			// Siniukow
	_types[2]  = "teposevent";	// EPOS
	
	_chain = new TChain();

	_therminatorParticle = 0;
}

TChainProxy::~TChainProxy()
{
	// delete _chain;
}

Int_t TChainProxy::Add(const char *fileName)
{
	// event type recognition on first file
	if(_type == NONE)
	{
		TFile *file;
		TTree 
			*therminatorTree = 0,
			*siniukowTree = 0,
			*eposTree = 0;

		file = new TFile(fileName);
		file->GetObject(_types[0], therminatorTree);
		file->GetObject(_types[1], siniukowTree);
		file->GetObject(_types[2], eposTree);

		if(therminatorTree)
			_type = Therminator;
		else if(siniukowTree)
			_type = Siniukow;
		else if(eposTree)
			_type = EPOS;
		else
		{
			std::cout << "Unknown event file: <" << fileName << "> !!!" << std::endl;
			delete file;
			return 0;
		}
		_chain->SetName(_types[_type]);
		delete file;
	}

	if(_type == EPOS)
	{// loading EPOS files
		if(_EPOSeventIndices.size() == 0) // loading only once
			_loadEPOSfilesChain();
		return 1;
	}
	else
		return _chain->Add(fileName);
}
// WARNING!
Int_t TChainProxy::GetEntry(Long64_t entryNumber)
{
	unsigned long  int previouslyLoadedParticleCount = 0, index = 0;
	long int iMax = _iEvent - 1 ;

	// if Therminator just simply return the number
	if(_type == Therminator)
		return _chain->GetEntry(entryNumber);
	// count already read
	for(int i = 0; i < iMax; ++i)
	{
		previouslyLoadedParticleCount += _particlesPerEvent[i];
		// std::cout << "[l00p] \t" << i << "\tval:\t" << _particlesPerEvent[i] << std::endl;
	}

	// entryNumber to index in _particlesPerEvent vector conversion
	index = entryNumber - previouslyLoadedParticleCount;

	if(index < 0 || index > _particlesPerEvent.size())
	{
		std::cerr << "CRITICAL ERROR! INDEX OUT OF BOUNDS!" << std::endl
		<< "\tINDEX:\t" << index << "\t_particlesPerEvent.size():\t" << _particlesPerEvent.size() << std::endl;
	}

	// std::cout << "[GetEntry] entry:\t" << entryNumber << " evt:\t" << _iEvent << " p:\t" << previouslyLoadedParticleCount << " i:\t" << index 
	// << " vector size:\t" << _EPOSParticles.size() << std::endl;

	if(_type == Siniukow)
	{
		if(index >= _siniukowParticles.size())
			_isSiniukowLoaded = false;
		_loadSiniukow();
		*_therminatorParticle = _siniukowParticles[index];
	}
	else if(_type == EPOS)
	{
		if(index >= _EPOSParticles.size())
			_isEPOSLoaded = false;
		_loadEPOS();
		*_therminatorParticle = _EPOSParticles[index];
	}

	return 1;
}

Long64_t TChainProxy::GetEntries()
{
	unsigned int i, j, counter, numberOfEvents;
	Int_t particlesCount;
	Long64_t totalParticlesCount = 0;

	int pionCounter = 0;
	Int_t id[_maxParticleCount];
	
	Int_t eposStatus[_maxParticleCount]; // required for particle selection

	if(_type == Therminator)
		return _chain->GetEntries();
	else
	{ // counting all particles in all events
		_chain->SetBranchStatus("*",1);
		if(_type == Siniukow)
			_chain->SetBranchAddress("npart", &particlesCount);
		else if(_type == EPOS)
		{
			_chain->SetBranchAddress("np", &particlesCount);
			_chain->SetBranchAddress("ist", eposStatus);
		}

		if(_type == EPOS)
			numberOfEvents = _EPOSeventIndices.size();
		else
			numberOfEvents = _chain->GetEntries();

		for(i = 0; i < numberOfEvents; ++i)
		{
			if(_type == EPOS)
				_chain->GetEntry(_EPOSeventIndices[i]);
			else
				_chain->GetEntry(i);

			if(_type == Siniukow)
			{
				_particlesPerEvent.push_back(particlesCount); // all particles are taken into account
				totalParticlesCount += particlesCount;
			}
			else if(_type == EPOS)
			{
				counter = 0;
				for(j = 0; j < particlesCount; ++j)
					if(eposStatus[j] == 0) // only last generation hadrons
					{
						++counter;
					}
				_particlesPerEvent.push_back(counter);
				totalParticlesCount += counter;
			}
			// std::cout << "[GetEntries] e:\t" << i << " p:\t" << _particlesPerEvent[i] << std::endl;
		}
	}
	if(_type == EPOS)
		std::cout << "EPOS: entries: " << totalParticlesCount << " @ evts: " << numberOfEvents << " Last chainevent: " << _chain->GetEntries() << std::endl;
	return totalParticlesCount;
}

Int_t TChainProxy::SetBranchAddress(const char *ignored, ParticleCoor *particle)
{
	_therminatorParticle = particle;
	if(_type == Therminator)
		_chain->SetBranchAddress("particle", _therminatorParticle);
	else if (_type == Siniukow)
		;	// nothing here
	else if (_type == EPOS)
		; 	// nothing here
	else
	{
		std::cout << "Unknown event file!!!" << std::endl;
		return 1;
	}
	return 1;
}

void TChainProxy::_loadSiniukow()
{
	SiniukowEvent event;
	ParticleCoor particle;
	unsigned int iParticle, eid;

	if(_type != Siniukow)
		return; // ONLY for Siniukow
	if(_isSiniukowLoaded)
		return; // event files are already loaded

	_chain->SetBranchStatus("*",1);
	_chain->SetBranchAddress("npart",&event.npart);
	_chain->SetBranchAddress("id", event.id);
	_chain->SetBranchAddress("mid", event.mid );
	_chain->SetBranchAddress("x", event.x);
	_chain->SetBranchAddress("y", event.y);
	_chain->SetBranchAddress("z", event.z);
	_chain->SetBranchAddress("t", event.t);
	_chain->SetBranchAddress("px", event.px);
	_chain->SetBranchAddress("py", event.py);
	_chain->SetBranchAddress("pz", event.pz);
	_chain->SetBranchAddress("E", event.E);

	_totalEventCount = _chain->GetEntries();
	_siniukowParticles.clear();

	if(_iEvent < _totalEventCount)
	{
		eid = 0;
		_chain->GetEntry(_iEvent);
		for (iParticle = 0; iParticle < event.npart; ++iParticle)
		{
			particle.mass = TMath::Sqrt(
				TMath::Power(event.E[iParticle]/TMath::C(),2)
				- TMath::Power(event.px[iParticle],2)
				- TMath::Power(event.py[iParticle],2)
				- TMath::Power(event.pz[iParticle],2)
			);
			particle.t = event.t[iParticle];
			particle.x = event.x[iParticle];
			particle.y = event.y[iParticle];
			particle.z = event.z[iParticle];
			particle.e = event.E[iParticle];
			particle.px = event.px[iParticle];
			particle.py = event.py[iParticle];
			particle.pz = event.pz[iParticle];
			particle.decayed = 0;
			particle.pid = event.id[iParticle];
			particle.fatherpid = event.mid[iParticle];
			particle.rootpid = 0;
			particle.eid = ++eid;
			particle.fathereid = 0;
			particle.eventid = _iEvent;

			_siniukowParticles.push_back(particle);
		}
	}
	_isSiniukowLoaded = true;
	++_iEvent;
}

void TChainProxy::_loadEPOS()
{
	EPOSEvent event;
	ParticleCoor particle;
	unsigned int iParticle, eid;

	if(_type != EPOS)
		return; // ONLY for EPOS
	if(_isEPOSLoaded)
		return; // event files are already loaded

	_chain->SetBranchStatus("*",1);
	_chain->SetBranchAddress("np", &event.np);
	_chain->SetBranchAddress("bim", &event.bim);
	_chain->SetBranchAddress("zus", event.zus);
	_chain->SetBranchAddress("px", event.px);
	_chain->SetBranchAddress("py", event.py);
	_chain->SetBranchAddress("pz", event.pz);
	_chain->SetBranchAddress("e", event.e);
	_chain->SetBranchAddress("x", event.x);
	_chain->SetBranchAddress("y", event.y);
	_chain->SetBranchAddress("z", event.z);
	_chain->SetBranchAddress("t", event.t);
	_chain->SetBranchAddress("id", event.id);
	_chain->SetBranchAddress("ist", event.ist);
	_chain->SetBranchAddress("ity", event.ity);
	_chain->SetBranchAddress("ior", event.ior);
	_chain->SetBranchAddress("jor", event.jor);

	_totalEventCount = _EPOSeventIndices.size();
	_EPOSParticles.clear();

	// std::cout << "Loading EPOS Event nr " << _iEvent << " / " << _totalEventCount << std::endl;

	if(_iEvent < _totalEventCount)
	{
		eid = 0;
		if(_iEvent >= _EPOSeventIndices.size())
		{
			std::cout << "[TChainProxy] Index " << _iEvent << " is out of range: " << _EPOSeventIndices.size() << std::endl;
			return;
		}
		if(_chain->GetEntry(_EPOSeventIndices[_iEvent]) == 0)
		{
			std::cout << "[TChainProxy] Failed to load event: " << _iEvent << " - eposEvent: " << _EPOSeventIndices[_iEvent] << std::endl;
			return;
		}
		for (iParticle = 0; iParticle < event.np; ++iParticle)
		{
			if(event.ist[iParticle] != 0) // only last generation hadrons
				continue;

			particle.mass = TMath::Sqrt(
				TMath::Power(event.e[iParticle]/TMath::C(),2)
				- TMath::Power(event.px[iParticle],2)
				- TMath::Power(event.py[iParticle],2)
				- TMath::Power(event.pz[iParticle],2)
			);
			particle.t = event.t[iParticle];
			particle.x = event.x[iParticle];
			particle.y = event.y[iParticle];
			particle.z = event.z[iParticle];
			particle.e = event.e[iParticle];
			particle.px = event.px[iParticle];
			particle.py = event.py[iParticle];
			particle.pz = event.pz[iParticle];
			particle.decayed = 0;
			particle.pid = _translateEPOSpidToTherminator(event.id[iParticle]);
			if(event.ior[iParticle] == 0 && event.jor[iParticle] == 0)
			{
				particle.fatherpid = 0;
				particle.rootpid = 0;
				particle.fathereid = 0;
			}
			else
			{	particle.fatherpid = _translateEPOSpidToTherminator(event.id[event.ior[iParticle]]);
				particle.rootpid = 0;
				particle.fathereid = _translateEPOSpidToTherminator(event.ior[iParticle]);
			}
			particle.eid = ++eid;
			particle.eventid = _iEvent;

			_EPOSParticles.push_back(particle);
		}
	}
	_isEPOSLoaded = true;
	++_iEvent;
}

void TChainProxy::_loadEPOSfilesChain()
{
	char buffer[256];
	unsigned int eventNum, filesToLoad, i,j;
	std::ifstream eventNumbers("events.centrality", std::ifstream::in);
	std::ifstream eventFiles("eventFiles.list", std::ifstream::in);

	if(!eventNumbers.is_open())
	{
		std::cout << "Could not open \"events.centrality\" file!!! Quitting." << std::endl;
		exit(1);
	}

	if(!eventFiles.is_open())
	{
		std::cout << "Could not open \"eventFiles.list\" file!!! Quitting." << std::endl;
		exit(1);
	}

	j = 0;
	while(eventNumbers.good())
	{
		for(i=0; i < 256; ++i)
			buffer[i] = '\0';
    	eventNumbers >> buffer;
    	if(buffer[0] == '\0')
    		continue;
    	std::stringstream(buffer) >> eventNum;
    	_EPOSeventIndices.push_back(eventNum);
    	// std::cout << ++j << " " << eventNum << std::endl;
    }

    filesToLoad = _EPOSeventIndices.back() / _EPOSeventsPerFile + 1;


	for(i =0; (i < filesToLoad) && eventFiles.good(); ++i)
	{
		for(j=0; j < 256; ++j)
			buffer[j] = '\0';
    	eventFiles >> buffer;
    	if(buffer[0] == '\0')
    		continue;
//    	std::cout << "[TChainProxy|EPOS] Loading:\t\"" << buffer << "\"" << std::endl;
    	if(_chain->Add(buffer) == 0)
	    	std::cout << "[TChainProxy|EPOS]  Problem loading file: \t\"" << buffer << "\" !!!" << std::endl;
    }
}

Int_t TChainProxy::_translateEPOSpidToTherminator(Int_t pid)
{
	Bool_t isAntiParticle = false;
	if(pid < 0 )
	{	
		isAntiParticle = true;
		pid = -pid;
	}

	switch(pid)
	{
		case 10:
			pid = 22;
			break;
		case 110:
			pid = 111;
			break;
		case 120:
			pid = 211;
			break;
		case -120:
			pid = -211;
			break;
		case 220:
			pid = 221;
			break;
		case -220:
			pid = -221;
			break;
		case 130:
			pid = 321;
		case -130:
			pid = -321;
			break;
		case 230:
			pid = 311;
			break;
		case -230:
			pid = -311;
			break;
		case 330:
			pid = 331;
			break;
		case 140:
			pid = 1234;
			break;
		case -140:
			pid = 1233;
			break; 
		case 240:
			pid = 1232;
			break; 
		case -240:
			pid = 1231;
			break; 
		// case 340:
		// 	pid = 
		// 	break; 
		// case -340:
		// 	pid = 
		// 	break; 
		case 111:
			pid = 113;
			break;
		case 121:
			pid = 213;
			break;
		case -121:
			pid = -213;
			break;
		case 221:
			pid = 223;
			break;
		case 131:
			pid = 323;
			break;
		case -131:
			pid = -323;
			break;
		case 231:
			pid = 313;
			break;
		case -231:
			pid = -313;
			break;
		case 331:
			pid = 333;
			break;
		case 141:
			pid = 1234;
			break;
		case -141:
			pid = 1233;
			break;
		case 241:
			pid = 1232;
			break;
		case -241:
			pid = 1231;
			break;
		case 441:
			pid = 443;
			break;
		case 1120:
			pid = 2212;
			break;
		case 1220:
			pid = 2112;
			break;
		case 1130:
			pid = 3222;
			break;
		case 1230:
			pid = 3212;
			break;
		case 2130:
			pid = 3122;
			break;
		case 2230:
			pid = 3112;
			break;
		case 1330:
			pid = 3322;
			break;
		case 2330:
			pid = 3312;
			break;
		// 1140 - 3440 not found
		case 1111:
			pid = 11114;
			break;
		case 1121:
			pid = 12114;
			break;
		case 1221:
			pid = 12214;
			break;
		case 2221:
			pid = 12224;
			break;
		case 1131:
			pid = 3224;
			break;
		case 1231:
			pid = 3214;
			break;
		case 2231:
			pid = 3114;
			break;
		case 1331:
			pid = 3314;
			break;
		case 2331:
			pid = 3324;
			break;
		case 3331:
			pid = 227;
			break;
		// 1141 - 4441 - not found
		case 332:
			pid = 9010221;
			break;
		case 112:
			pid = 9000111;
			break;
		case 122:
			pid = 9000211;
			break;
		case 1112:
			pid = 1112;
			break;
		case 1113:
			pid = 11114;
			break;
		case 1114:
			pid = 11116;
			break;
		case 2222:
			pid = 2222;
			break;
		case 2223:
			pid = 12224;
			break;
		case 2224:
			pid = 22224;
			break;
		case 1122:
			pid = 12212;
			break;
		case 1123:
			pid = 22212;
			break;
		case 1124:
			pid = 1212;
			break;
		case 1125:
			pid = 32112;
			break;
		case 1126:
			pid = 12114;
			break;
		case 1127:
			pid = 42112;
			break;
		case 1128:
			pid = 11216;
			break;
		case 1222:
			pid = 12112;
			break;
		case 1223:
			pid = 22122;
			break;
		case 1224:
			pid = 2122;
			break;
		case 1225:
			pid = 32212;
			break;
		case 1226:
			pid = 12214;
			break;
		case 1227:
			pid = 42212;
			break;
		case 1228:
			pid = 22214;
			break;
		case 1233:
			pid = 13122;
			break;
		case 1234:
			pid = 3124;
			break;
		case 1236:
			pid = 13224;
			break;
		case 1237:
			pid = 3226;
			break;
		case 1238:
			pid = 13126;
			break;
		case 1239:
			pid = 23224;
			break;
		case 1132:
			pid = 13114;
			break;
		case 1133:
			pid = 3116;
			break;
		case 1134:
			pid = 23114;
			break;
		case 2232:
			pid = 13214;
			break;
		case 2233:
			pid = 3216;
			break;
		case 2234:
			pid = 23214;
			break;
	}
	if(isAntiParticle)
		pid = (-1)*pid;
	return pid;
}
