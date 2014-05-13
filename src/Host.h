#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include <list>
#include "EventQueue.hpp"
#include "Strain.h"
#include "Gene.h"
#include "zppsim_util.hpp"
#include "SimParameters.h"
#include "Infection.h"
#include "ImmuneHistory.h"
#include "Database.hpp"

#define WAITING_STAGE (std::numeric_limits<size_t>::max())

class Host;
class Population;
class Infection;

class DeathEvent : public zppsim::Event
{
public:
	DeathEvent(Host * hostPtr);
	virtual void performEvent(zppsim::EventQueue & queue);
	
	Host * hostPtr;
};


class Host
{
friend class Population;
friend class DeathEvent;
friend class Infection;
friend class Immunity;
public:
	size_t const id;
	
	Host(Population * popPtr, size_t id, double birthTime, double deathTime, zppdata::DBTable * table);
	
	void prepareToDie();
	
	void transmitTo(Host & dstHost);
	
	void receiveInfection(StrainPtr & strain);
	
	void updateInfectionRates();
	
	double getAge();
	
	size_t getActiveInfectionCount();
	std::vector<GenePtr> getActiveInfectionGenes();
	std::vector<size_t> getActiveInfectionGeneIds();
	
	size_t getActiveInfectionImmunityCount();
	size_t getActiveInfectionClinicalImmunityCount();
	
	void clearInfection(std::list<Infection>::iterator infectionItr);
	
	double getTime();
	zppsim::rng_t * getRngPtr();
	
	SimParameters * getSimulationParametersPtr();
	PopulationParameters * getPopulationParametersPtr();
	
	void addEvent(zppsim::Event * event);
	void removeEvent(zppsim::Event * event);
	void setEventRate(zppsim::RateEvent * event, double rate);
	
	void writeInfections(DBTable * table);
	
	std::string toString();
private:
	Population * popPtr;
	double const birthTime;
	double const deathTime;
	
	size_t nextInfectionId;
	
	// Linked list of current infections
	std::list<Infection> infections;
	
	std::unique_ptr<DeathEvent> deathEvent;
	
	// Two sets of immune history (regular & "clinical")
	ImmuneHistory immunity;
	ImmuneHistory clinicalImmunity;
};

#endif
