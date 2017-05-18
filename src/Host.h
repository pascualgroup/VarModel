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
#include "zppdb.hpp"
#include "DatabaseTypes.h"

#define WAITING_STAGE (std::numeric_limits<int64_t>::max())

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
friend class Simulation;
friend class Population;
friend class DeathEvent;
friend class Infection;
friend class ImmuneHistory;
public:
	int64_t const id;
	
	Host(
		Population * popPtr, int64_t id, double birthTime, double deathTime,
		bool writeToDatabase,
		Database & db,
		zppdb::Table<HostRow> & table
	);
	
    //whether to track this hosts' entire infection histroy
    bool toTrack = false;
    
	void prepareToDie();
    
    double moiRegulate(Host & dstHost);
	
	void transmitTo(Host & dstHost);

    void transmitMSTo(Host & dstHost);

	void receiveInfection(StrainPtr & strain);
	
    void receiveInfection(StrainPtr & strain, GenePtr & msPtr);
	
    void updateInfectionRates();
	
	double getAge();
	
	int64_t getActiveInfectionCount();
	std::vector<GenePtr> getActiveInfectionGenes();
	std::vector<int64_t> getActiveInfectionGeneIds();
    
    void MDAClearInfection();
	
	int64_t getActiveInfectionImmunityCount();
	
    void gainAlleleImmunity(GenePtr genePtr);
	void clearInfection(std::list<Infection>::iterator infectionItr);
    void hstMutateStrain(std::list<Infection>::iterator infectionItr);
    void RecombineStrain(std::list<Infection>::iterator infectionItr);
    
    void getSelectionMode(GenePtr genePtr, bool clearInfection);
	void microsatMutate(std::list<Infection>::iterator infectionItr);
    
	double getTime();
	zppsim::rng_t * getRngPtr();
	  
	SimParameters * getSimulationParametersPtr();
	PopulationParameters * getPopulationParametersPtr();
	
	void addEvent(zppsim::Event * event);
	void removeEvent(zppsim::Event * event);
	void setEventRate(zppsim::RateEvent * event, double rate);
	
	void writeInfections(Database & db, Table<InfectionRow> & table, Table<StrainRow> & strainsTable,Table<GeneRow> & GeneTable,Table<LociRow> & LociTable);
	void writeInfections(int64_t transmissionId, Database & db, Table<TransmissionInfectionRow> & table);
	
	std::string toString();
    
    void writeToCheckpoint(Database & cpdb);
private:
	Population * popPtr;
	double const birthTime;
	double const deathTime;
	
	int64_t nextInfectionId;
	
    double MDAEndTime = 0;
    
	// Linked list of current infections
	std::list<Infection> infections;
	
	std::unique_ptr<DeathEvent> deathEvent;
	
	ImmuneHistory immunity;
};

#endif
