#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;
using namespace zppsim;

Host::Host(Population * popPtr, size_t id, double birthTime, double deathTime) :
	popPtr(popPtr), id(id),
	birthTime(birthTime), deathTime(deathTime), nextInfectionId(0),
	deathEvent(new DeathEvent(this))
{
//	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	addEvent(deathEvent.get());
}

double Host::getAge()
{
	return getTime() - birthTime;
}

size_t Host::infectionCount()
{
	return infections.size();
}

bool Host::isImmune(GenePtr gene)
{
	return immunity.find(gene) != immunity.end();
}

void Host::prepareToDie()
{
//	cerr << popPtr->getTime() << ", host going to die: " << toString() << '\n';
	
	// Remove all events
	for(auto & infection : infections) {
		infection.prepareToEnd();
	}
	for(auto itr = immunityLossEvents.begin(); itr != immunityLossEvents.end(); itr++) {
		removeEvent(itr->second.get());
	}
	removeEvent(deathEvent.get());
}

void Host::transmitTo(Host & dstHost)
{
	rng_t * rngPtr = getRngPtr();
	
	if(infections.size() == 0) {
	//	cerr << "No infections to transmit" << endl;
		return;
	}
//	cerr << "Transmitting to " << dstHost.popPtr->id << ", " << dstHost.id << '\n';
	
	// Get some current infections according to transmission probability
	vector<StrainPtr> originalStrains;
	for(auto infItr = infections.begin(); infItr != infections.end(); infItr++) {
		if(infItr->isActive()) {
			bernoulli_distribution flipCoin(infItr->transmissionProbability());
			if(flipCoin(*rngPtr)) {
				originalStrains.push_back(
					popPtr->simPtr->mutateStrain(infItr->strainPtr)
				);
			}
		}
	}
	
	vector<StrainPtr> strainsToTransmit;
	if(originalStrains.size() > 1) {
		// Take each original strain with probability 1.0 - pRecombinant
		double pRecombinant = popPtr->simPtr->parPtr->pRecombinant;
		assert(pRecombinant >= 0.0 && pRecombinant <= 1.0);
		strainsToTransmit = drawMultipleBernoulliRandomSubset(
			*rngPtr, originalStrains, 1.0 - pRecombinant, false
		);
		
		// Complete a set of size originalStrains.size() using recombinants
		size_t nRecombinants = originalStrains.size() - strainsToTransmit.size();
		for(size_t i = 0; i < nRecombinants; i++) {
			uniform_int_distribution<size_t> indDist(0, originalStrains.size() - 1);
			size_t ind1 = indDist(*rngPtr);
			size_t ind2 = indDist(*rngPtr);
			strainsToTransmit.push_back(
				popPtr->simPtr->recombineStrains(
					originalStrains[ind1],
					originalStrains[ind2]
				)
			);
		}
	}
	else {
		strainsToTransmit = originalStrains;
	}
	
	// Transmit all strains
	for(auto & strainPtr : strainsToTransmit) {
		dstHost.receiveInfection(strainPtr);
	}
}

void Host::receiveInfection(StrainPtr & strainPtr)
{
	assert(strainPtr->size() > 0);
	
	popPtr->countTransmission();
	
	rng_t * rngPtr = getRngPtr();
	double time = popPtr->getTime();
	
	double tLiverStage = popPtr->simPtr->parPtr->tLiverStage;
	size_t initialGeneIndex = tLiverStage == 0 ? 0 : LIVER_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, initialGeneIndex);
	
	
	// Get an iterator to the infection so that events can
	// point back to their infection
	list<Infection>::reverse_iterator reverseItr = infections.rbegin();
	reverseItr++;
	list<Infection>::iterator infectionItr = reverseItr.base();
	
	// If starting in liver stage, create a fixed-time transition event
	// (liver stage -> first gene not yet active)
	if(initialGeneIndex == LIVER_STAGE) {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(infectionItr, time + tLiverStage)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	// Otherwise create a rate-based transition event
	// (first gene not yet active -> first gene active)
	else {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(
				infectionItr,
				infectionItr->transitionRate(),
				time,
				*rngPtr
			)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	
	// Create a clearance event
	// (rate will depend on state as determined in clearanceRate()
	// and may be zero)
	infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
		new ClearanceEvent(
			infectionItr, infectionItr->clearanceRate(), time, *rngPtr
		)
	);
	addEvent(infectionItr->clearanceEvent.get());
	
//	cerr << time << ": " << infectionItr->toString() << " begun" << '\n';
}

void Host::gainImmunity(GenePtr genePtr)
{
	if(immunity.find(genePtr) == immunity.end()) {
		immunity.insert(genePtr);
		
		assert(immunityLossEvents.find(genePtr) == immunityLossEvents.end());
		ImmunityLossEvent * ilEvent = new ImmunityLossEvent(
			this, genePtr, genePtr->immunityLossRate, getTime()
		);
		immunityLossEvents[genePtr] = unique_ptr<ImmunityLossEvent>(ilEvent);
		addEvent(ilEvent);
		
		updateInfectionRates();
	
	//	cerr << getTime() << ": " << toString() << " gained immunity to " << genePtr->toString() << '\n';
	}
	else {
	//	cerr << getTime() << ": " << toString() << " already had immunity to " << genePtr->toString() << '\n';
	}
}

void Host::loseImmunity(GenePtr genePtr)
{
	// Remove immunity
	auto itr1 = immunity.find(genePtr);
	assert(itr1 != immunity.end());
	immunity.erase(itr1);
	
	// Remove immunity loss event
	auto itr2 = immunityLossEvents.find(genePtr);
	assert(itr2 != immunityLossEvents.end());
	ImmunityLossEvent * ilEvent = itr2->second.get();
	removeEvent(ilEvent);
	immunityLossEvents.erase(itr2);
	
	updateInfectionRates();
	
//	cerr << getTime() << ": " << toString() << " lost immunity to " << genePtr->toString() << '\n';
}

void Host::updateInfectionRates()
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		if(itr->geneIndex != LIVER_STAGE) {
			itr->updateClearanceRate();
			itr->updateTransitionRate();
		}
	}
}

void Host::clearInfection(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	double time = popPtr->getTime();
//	cerr << time << ": " << infectionItr->toString() << " clearing" << '\n';
	
	// Gain immunity to active gene
	if(infectionItr->active) {
		gainImmunity(infectionItr->getCurrentGene());
	}
	
	// Remove infection
	infectionItr->prepareToEnd();
	infections.erase(infectionItr);
}

double Host::getTime()
{
	return popPtr->getTime();
}

rng_t * Host::getRngPtr()
{
	return popPtr->rngPtr;
}

SimParameters * Host::getSimulationParametersPtr()
{
	return popPtr->simPtr->parPtr;
}

PopulationParameters * Host::getPopulationParametersPtr()
{
	return popPtr->parPtr;
}

void Host::addEvent(Event * event)
{
	popPtr->addEvent(event);
}

void Host::removeEvent(Event * event)
{
	popPtr->removeEvent(event);
}

void Host::setEventRate(RateEvent * event, double rate)
{
	popPtr->setEventRate(event, rate);
}

std::string Host::toString()
{
	return popPtr->toString() + ".h" + strprintf("%u", id);
}

Infection::Infection(Host * hostPtr, size_t id, StrainPtr & strainPtr, size_t initialGeneIndex) :
	hostPtr(hostPtr), id(id), strainPtr(strainPtr),
	geneIndex(initialGeneIndex), active(false)
{
}

void Infection::prepareToEnd()
{
	hostPtr->removeEvent(transitionEvent.get());
	hostPtr->removeEvent(clearanceEvent.get());
}

GenePtr Infection::getCurrentGene()
{
	assert(geneIndex != LIVER_STAGE);
	return strainPtr->getGene(geneIndex);
}

void Infection::performTransition()
{
	if(geneIndex == LIVER_STAGE) {
		assert(!active);
		geneIndex = 0;
	}
	else if(active) {
		assert(geneIndex != strainPtr->size() - 1);
		hostPtr->gainImmunity(strainPtr->getGene(geneIndex));
		geneIndex++;
		active = false;
	}
	else {
		active = true;
	}
	
	updateTransitionRate();
	updateClearanceRate();
	
	double time = hostPtr->getTime();
//	cerr << time << ": " << toString() << " transitioned to " << geneIndex << "(" << getCurrentGene()->toString() << "), " << (active ? "active" : "not yet active") << '\n';
}

void Infection::updateTransitionRate()
{
	hostPtr->setEventRate(transitionEvent.get(), transitionRate());
}

double Infection::transitionRate()
{
	assert(geneIndex != LIVER_STAGE);
	
	if(!active) {
		return activationRate();
	}
	else {
		return deactivationRate();
	}
}

void Infection::updateClearanceRate()
{
	hostPtr->setEventRate(clearanceEvent.get(), clearanceRate());
}

double Infection::activationRate()
{
//	GenePtr gene = getCurrentGene();
//	bool immune = hostPtr->isImmune(gene);
//	double age = hostPtr->getAge();
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	return simParPtr->withinHost.activationRate;
}

double Infection::deactivationRate()
{
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	GenePtr gene = getCurrentGene();
	bool immune = hostPtr->isImmune(gene);
	if(immune) {
		return simParPtr->withinHost.deactivationRateImmune;
	}
	else {
		return simParPtr->withinHost.deactivationRateNotImmune;
	}
}

double Infection::clearanceRate()
{
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	// Liver stage
	if(geneIndex == LIVER_STAGE) {
		return 0.0;
	}
	// Gene active
	else if(!active) {
		// Before first activation
		if(geneIndex == 0) {
			return simParPtr->withinHost.clearanceRateInitial;
		}
		// Between activations
		else {
			return simParPtr->withinHost.clearanceRateMidCourse;
		}
	}
	// Gene active
	else {
		return 0.0;
	}
}

bool Infection::isActive()
{
	return active;
}

double Infection::transmissionProbability()
{
	assert(active);
	
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	GenePtr genePtr = getCurrentGene();
	
	double p = genePtr->transmissibility;
	assert(p > 0.0 && p < 1.0);
	
	// TODO: reduce by active infections or all infections?
	if(simParPtr->transmission.coinfectionReducesTransmission) {
		p /= hostPtr->infectionCount();
	}
	return p;
}

std::string Infection::toString()
{
	return hostPtr->toString() + ".i" + strprintf("%u", id);
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->prepareToDie();
	hostPtr->popPtr->removeHost(hostPtr);
}

ImmunityLossEvent::ImmunityLossEvent(Host * hostPtr, GenePtr genePtr,
	double rate, double initTime) :
	RateEvent(rate, initTime, *hostPtr->getRngPtr()),
	hostPtr(hostPtr), genePtr(genePtr)
{
}

void ImmunityLossEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->loseImmunity(genePtr);
}

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	RateEvent(time), infectionItr(infectionItr)
{
}

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	RateEvent(rate, time, rng), infectionItr(infectionItr)
{
}

TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	InfectionProcessEvent(infectionItr, time)
{
}


TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

ClearanceEvent::ClearanceEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

void TransitionEvent::performEvent(zppsim::EventQueue &queue)
{
	// If this is the final deactivation, then it's equivalent to clearing
	if(infectionItr->active
		&& infectionItr->geneIndex == infectionItr->strainPtr->size() - 1
	) {
		infectionItr->hostPtr->clearInfection(infectionItr);
	}
	// Otherwise, actually perform a transition
	else {
		infectionItr->performTransition();
	}
}

void ClearanceEvent::performEvent(zppsim::EventQueue &queue)
{
	infectionItr->hostPtr->clearInfection(infectionItr);
}
