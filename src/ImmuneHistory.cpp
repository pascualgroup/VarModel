//
//  ImmuneHistory.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

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


/*** ImmuneHistory function implementations ***/

ImmuneHistory::ImmuneHistory(Host * hostPtr, bool clinical) : hostPtr(hostPtr), clinical(clinical)
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
