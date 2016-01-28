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
/*
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
*/
/*** ImmuneHistory function implementations ***/

ImmuneHistory::ImmuneHistory(Host * hostPtr, bool clinical, int64_t const locusNumber, double infectionTimesToImmune) : hostPtr(hostPtr), clinical(clinical),locusNumber(locusNumber), infectionTimesToImmune(infectionTimesToImmune)
{
}

void ImmuneHistory::gainImmunity(GenePtr genePtr)
{
	if(genes.find(genePtr) == genes.end()) {
		genes.insert(genePtr);
		
		assert(lossEvents.find(genePtr) == lossEvents.end());
		double immunityLossRate = clinical ?
			genePtr->clinicalImmunityLossRate :
			genePtr->immunityLossRate;
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

double ImmuneHistory::checkGeneralImmunity() {
    //can add this as a parameter later
    if (infectedTimes >= infectionTimesToImmune) {
        return 1.0;
    }else{
            return infectedTimes/infectionTimesToImmune;
    }
}

void ImmuneHistory::gainAlleleImmunity(GenePtr genePtr,bool writeToDatabase,Database & db,zppdb::Table<AlleleImmunityRow> & table)
{
    vector<int64_t> geneAlleles = genePtr->Alleles;
    if(immuneAlleles.empty()) { //construct a list with immuned alleles
        for (int64_t i=0; i<locusNumber;i++) {
            std::unordered_map<int64_t,int64_t> tempMap({{geneAlleles[i],1}});
            immuneAlleles.push_back(tempMap);
        }
    }else{
        for (int64_t i=0; i<locusNumber;i++) {
            if(immuneAlleles[i].find(geneAlleles[i])==immuneAlleles[i].end()) {
                immuneAlleles[i][geneAlleles[i]] = 1;
                if(writeToDatabase) {
                    AlleleImmunityRow row;
                    row.time = hostPtr->getTime();
                    row.hostId = hostPtr->id;
                    row.locusIndex = i;
                    row.alleleId = geneAlleles[i];
                    db.insert(table, row);
                }
                    
            }else{
                immuneAlleles[i][geneAlleles[i]] ++;
            }
        }
    }
    hostPtr->updateInfectionRates();
    
}

/*
void ImmuneHistory::setAlleleLossEvent(int64_t locusId, int64_t AlleleId,double lossrate) {
    assert(alleleLossEvents[locusId].find(AlleleId) == alleleLossEvents[locusId].end());
    AlleleImmuneLossEvent * ilEvent = new AlleleImmuneLossEvent(
                                                        this, locusId, AlleleId,lossrate, hostPtr->getTime());
    alleleLossEvents[locusId][AlleleId] = unique_ptr<AlleleImmuneLossEvent>(ilEvent);
    hostPtr->addEvent(ilEvent);
    
}
 */

double ImmuneHistory::checkGeneImmunity(GenePtr genePtr) {
    vector<int64_t> geneAlleles = genePtr->Alleles;
    int64_t immuneLevel = locusNumber;
    if(immuneAlleles.empty()) {
        return 0;
    }else{
        for (int64_t i=0; i<locusNumber;i++) {
            if(immuneAlleles[i].find(geneAlleles[i])==immuneAlleles[i].end()) {
                immuneLevel--;
            }
        }
        double immuneFraction = (double)immuneLevel/(double)locusNumber;
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

/*
void ImmuneHistory::loseAlleleImmune(int64_t & locusId,int64_t & AlleleId)
{
	// Remove allele immunity
	auto itr1 = immuneAlleles[locusId].find(AlleleId);
	assert(itr1 != immuneAlleles[locusId].end());
	immuneAlleles[locusId].erase(itr1);
	
	// Remove immunity loss event
	auto itr2 = alleleLossEvents[locusId].find(AlleleId);
	assert(itr2 != alleleLossEvents[locusId].end());
	AlleleImmuneLossEvent * ilEvent = itr2->second.get();
	hostPtr->removeEvent(ilEvent);
	alleleLossEvents[locusId].erase(itr2);
	
	hostPtr->updateInfectionRates();
}
*/

bool ImmuneHistory::isImmune(GenePtr genePtr)
{
	return genes.find(genePtr) != genes.end();
}

void ImmuneHistory::prepareToDie()
{
	for(auto itr = lossEvents.begin(); itr != lossEvents.end(); itr++) {
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
