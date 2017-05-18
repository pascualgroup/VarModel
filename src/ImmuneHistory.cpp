//
//  ImmuneHistory.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

//  change to immunity based on epitope/alleles, hqx

#include "ImmuneHistory.h"
#include "Host.h"

using namespace std;
using namespace zppsim;


/*** ImmunityLossEvent function implementations ***/

ImmunityLossEvent::ImmunityLossEvent(ImmuneHistory * immHistPtr, GenePtr genePtr,
	double rate, double initTime) :
	RateEvent(rate, initTime, *immHistPtr->hostPtr->getRngPtr()),
	immHistPtr(immHistPtr), genePtr(genePtr)
{
}

void ImmunityLossEvent::performEvent(zppsim::EventQueue & queue)
{
	immHistPtr->loseImmunity(genePtr);
}

/*** AlleleImmuneLossEvent function implementations ***/

AlleleImmuneLossEvent::AlleleImmuneLossEvent(ImmuneHistory * immHistPtr, int64_t & locusId, int64_t & AlleleId,
                                     double rate, double initTime) :
RateEvent(rate, initTime, *immHistPtr->hostPtr->getRngPtr()),
immHistPtr(immHistPtr), locusId(locusId),AlleleId(AlleleId)
{
}

void AlleleImmuneLossEvent::performEvent(zppsim::EventQueue & queue)
{
	immHistPtr->loseAlleleImmune(locusId,AlleleId);
}

/*** ImmuneHistory function implementations ***/

ImmuneHistory::ImmuneHistory(Host * hostPtr, int64_t const locusNumber, double infectionTimesToImmune) : hostPtr(hostPtr), locusNumber(locusNumber), infectionTimesToImmune(infectionTimesToImmune)
{
}

void ImmuneHistory::gainImmunity(GenePtr genePtr)
{
	if(genes.find(genePtr) == genes.end()) {
		genes.insert(genePtr);
		
		assert(lossEvents.find(genePtr) == lossEvents.end());
		double immunityLossRate = genePtr->immunityLossRate;
		ImmunityLossEvent * ilEvent = new ImmunityLossEvent(
			this, genePtr, immunityLossRate, hostPtr->getTime()
		);
		lossEvents[genePtr] = unique_ptr<ImmunityLossEvent>(ilEvent);
		hostPtr->addEvent(ilEvent);
		
		hostPtr->updateInfectionRates();
	}
}

void ImmuneHistory::gainGeneralImmunity() {
        infectedTimes++;
        hostPtr->updateInfectionRates();
}

double ImmuneHistory::checkGeneralImmunity(double a, double b) {
    //can add this as a parameter later
    if (infectedTimes<infectionTimesToImmune) {
        double y = exp(a-b*double(infectedTimes));
        return 1/y;
    }else{
        return 1;
    }
}

double ImmuneHistory::checkGeneralImmunity(std::vector<double> params, double & immuneRate) {
    //can add this as a parameter later
    if (infectedTimes<infectionTimesToImmune) {
        double y = params[1]*exp(-params[2]*double(infectedTimes))/pow((params[3]*double(infectedTimes)+1.0),params[3])+params[0];
        //cout<<infectedTimes<<" "<<y<<endl;
        return 1/y;
    }else{
        return immuneRate;
    }
}

int64_t ImmuneHistory::geneLocToHash(int64_t locusId, int64_t AlleleId) {
    int64_t a = AlleleId*100;
    a += locusId;
    //cout<<"geneHashVal "<<a<<endl;
    return (a);
}


void ImmuneHistory::gainAlleleImmunity(GenePtr genePtr,bool writeToDatabase,Database & db,zppdb::Table<AlleleImmunityRow> & table)
{
    double immunityLossRate = genePtr->immunityLossRate;
    vector<int64_t> geneAlleles = genePtr->Alleles;
    if(immuneAlleles.empty()) { //construct a list with immuned alleles
        for (int64_t i=0; i<locusNumber;i++) {
            std::unordered_map<int64_t,int64_t> tempMap({{geneAlleles[i],1}});
            immuneAlleles.push_back(tempMap);
            //set Allele loss event for each allele
            //rate inverse proportional to infected times
            setAlleleLossEvent(i,geneAlleles[i],immunityLossRate);
            //cout<<"addImmunity at locus "<<i<<" of allele "<<geneAlleles[i]<<endl;
        }
    }else{
        for (int64_t i=0; i<locusNumber;i++) {
            if(immuneAlleles[i].find(geneAlleles[i])==immuneAlleles[i].end()) {
                immuneAlleles[i][geneAlleles[i]] = 1;
                //cout<<"addImmunity at locus "<<i<<" of allele "<<geneAlleles[i]<<endl;
                if(writeToDatabase) {
                    AlleleImmunityRow row;
                    row.time = hostPtr->getTime();
                    row.hostId = hostPtr->id;
                    row.locusIndex = i;
                    row.alleleId = geneAlleles[i];
                    db.insert(table, row);
                }
                setAlleleLossEvent(i,geneAlleles[i],immunityLossRate);
                
            }else{
                immuneAlleles[i][geneAlleles[i]] ++;
                updateAlleleLossRate(i, geneAlleles[i],immunityLossRate/immuneAlleles[i][geneAlleles[i]]);
            }
        }
    }
    //hostPtr->updateInfectionRates();
    
}

void ImmuneHistory::setAlleleLossEvent(int64_t locusId, int64_t AlleleId,double lossrate) {
    double geneLocHashValue = geneLocToHash(locusId, AlleleId);
    assert(alleleLossEvents.find(geneLocHashValue) == alleleLossEvents.end());
    AlleleImmuneLossEvent * ilEvent = new AlleleImmuneLossEvent(
                                                        this, locusId, AlleleId,lossrate, hostPtr->getTime());
    alleleLossEvents[geneLocHashValue] = unique_ptr<AlleleImmuneLossEvent>(ilEvent);
    hostPtr->addEvent(ilEvent);
    
}

void ImmuneHistory::updateAlleleLossRate(int64_t & locusId, int64_t & AlleleId,double newRate) {
    auto itr2 = alleleLossEvents.find(geneLocToHash(locusId,AlleleId));
    assert(itr2 != alleleLossEvents.end());
    hostPtr->setEventRate(itr2->second.get(), newRate);
}

double ImmuneHistory::checkGeneImmunity(GenePtr genePtr) {
    vector<int64_t> geneAlleles = genePtr->Alleles;
    double immuneLevel = 0;
    double immuneTime = 1.0;
    std::unordered_map<int64_t,int64_t>::iterator it;
    if(immuneAlleles.empty()) {
        return 0;
    }else{
        for (int64_t i=0; i<locusNumber;i++) {
            it = immuneAlleles[i].find(geneAlleles[i]);
            if(it != immuneAlleles[i].end()) {
                if (it->second<=immuneTime) {
                    immuneLevel += 1/immuneTime*(it->second);
                }else{
                    immuneLevel += 1;
                }
            }        }
        double immuneFraction = immuneLevel/(double)locusNumber;
        //cout<<immuneFraction<<endl;
        return immuneFraction;
    }
}

void ImmuneHistory::loseImmunity(GenePtr genePtr)
{
	// Remove immunity
	auto itr1 = genes.find(genePtr);
	assert(itr1 != genes.end());
	genes.erase(itr1);
	
	// Remove immunity loss event
	auto itr2 = lossEvents.find(genePtr);
	assert(itr2 != lossEvents.end());
	ImmunityLossEvent * ilEvent = itr2->second.get();
	hostPtr->removeEvent(ilEvent);
	lossEvents.erase(itr2);
	
	hostPtr->updateInfectionRates();
}


void ImmuneHistory::loseAlleleImmune(int64_t & locusId,int64_t & AlleleId)
{
    //cout<<locusId<<":"<<AlleleId<<endl;
	// Remove allele immunity
	auto itr1 = immuneAlleles[locusId].find(AlleleId);
	assert(itr1 != immuneAlleles[locusId].end());
	immuneAlleles[locusId].erase(itr1);
	
	// Remove immunity loss event
	auto itr2 = alleleLossEvents.find(geneLocToHash(locusId,AlleleId));
	assert(itr2 != alleleLossEvents.end());
	AlleleImmuneLossEvent * ilEvent = itr2->second.get();
	hostPtr->removeEvent(ilEvent);
	alleleLossEvents.erase(itr2);
	hostPtr->updateInfectionRates();
}


bool ImmuneHistory::isImmune(GenePtr genePtr)
{
	return genes.find(genePtr) != genes.end();
}

void ImmuneHistory::prepareToDie()
{
	for(auto itr = lossEvents.begin(); itr != lossEvents.end(); itr++) {
		hostPtr->removeEvent(itr->second.get());
	}
	for(auto itr = alleleLossEvents.begin(); itr != alleleLossEvents.end(); itr++) {
		hostPtr->removeEvent(itr->second.get());
	}
    
}

void ImmuneHistory::write(Database & db, Table<ImmunityRow> & table)
{
	ImmunityRow row;
	row.time = hostPtr->getTime();
	row.hostId = hostPtr->id;
	for(auto & genePtr : genes) {
		row.geneId = genePtr->id;
		row.lossRate = lossEvents[genePtr]->getRate();
		db.insert(table, row);
	}
}

void ImmuneHistory::write(int64_t transmissionId, Database & db, Table<TransmissionImmunityRow> & table)
{
	TransmissionImmunityRow row;
	row.transmissionId = transmissionId;
	row.hostId = hostPtr->id;
	for(auto & genePtr : genes) {
		row.geneId = genePtr->id;
		row.lossRate = lossEvents[genePtr]->getRate();
		db.insert(table, row);
	}
}
