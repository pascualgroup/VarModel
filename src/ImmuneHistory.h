//
//  ImmuneHistory.h
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__ImmuneHistory__
#define __malariamodel__ImmuneHistory__

#include "EventQueue.hpp"
#include "zppsim_util.hpp"
#include <unordered_map>
#include <unordered_set>
#include "Gene.h"

class Host;
class ImmuneHistory;

class ImmunityLossEvent : public zppsim::RateEvent
{
public:
	ImmunityLossEvent(ImmuneHistory * immHistPtr, GenePtr genePtr, double rate, double initTime);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	ImmuneHistory * immHistPtr;
	GenePtr genePtr;
};

class ImmuneHistory
{
friend class ImmunityLossEvent;
public:
	ImmuneHistory(Host * hostPtr, bool clinical);
	
	void gainImmunity(GenePtr genePtr);
	void loseImmunity(GenePtr genePtr);
	bool isImmune(GenePtr genePtr);
	
	void prepareToDie();
	
	void write(Database & db, Table<ImmunityRow> & table);
	void write(int64_t transmissionId, Database & db, Table<TransmissionImmunityRow> & table);
	
	std::unordered_set<GenePtr> genes;
	std::unordered_map<GenePtr, std::unique_ptr<ImmunityLossEvent>> lossEvents;
private:
	Host * hostPtr;
	bool clinical;
};

#endif /* defined(__malariamodel__ImmuneHistory__) */
