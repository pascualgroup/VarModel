#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__


#include "SimParameters.h"
#include "Host.h"

#include "Database.hpp"
#include "zppsim_random.hpp"

#include <boost/array.hpp>
#include "zppdata_util.hpp"
#include "EventQueue.hpp"

class Simulation;

class BitingEvent : public zppsim::RateEvent
{
public:
	BitingEvent(Simulation * simPtr, uint32_t popId, zppsim::rng_t & rng);
	virtual void performEvent(zppsim::EventQueue & queue);
	
private:
	Simulation * simPtr;
	SimParameters * p;
	
	uint32_t popId;
};

class IntroductionEvent : public zppsim::RateEvent
{
public:
	IntroductionEvent(Simulation * simPtr, uint32_t popId, zppsim::rng_t & rng);
	virtual void performEvent(zppsim::EventQueue & queue);
	
private:
	Simulation * simPtr;
	SimParameters * p;
	
	uint32_t popId;
};

class Simulation
{
friend class BitingEvent;
friend class Host;
public:
	Simulation(SimParameters & params, zppdata::Database & db);
	
	void run();
	void runUntil(double time);
	void runOneEvent();
	
	uint64_t getNextHostId();
	double drawHostLifetime();
	
	double totalBitingRate(uint32_t popId);
	double totalIntroductionRate(uint32_t popId);
	
	void bite(uint32_t popId);
	void introduce(uint32_t popId);
	
	void markForDeath(Host & host);
	
	bool verifyState();
private:
	SimParameters * p;
	zppdata::Database * dbPtr;
	zppsim::rng_t rng;
	
	uint64_t nextHostId;
	
	std::vector<std::vector<std::unique_ptr<Host>>> hosts;
	std::unordered_map<uint64_t, std::pair<size_t, size_t>> hostIndexMap;
	
	Host * hostToDie;
	
	std::unique_ptr<zppsim::EventQueue> queuePtr;
	
	std::vector<std::unique_ptr<BitingEvent>> bitingEvents;
	std::vector<std::unique_ptr<IntroductionEvent>> introductionEvents;
};

#endif /* defined(__malariamodel__Simulation__) */
