//
//  ImmuneHistory.h
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

//  change to immunity based on epitope/alleles, hqx

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


class AlleleImmuneLossEvent : public zppsim::RateEvent
{
public:
	AlleleImmuneLossEvent(ImmuneHistory * immHistPtr, int64_t & locusId, int64_t & AlleleId, double rate, double initTime);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	ImmuneHistory * immHistPtr;
	int64_t locusId;
    int64_t AlleleId;
};


class ImmuneHistory
{
friend class ImmunityLossEvent;
friend class AlleleImmuneLossEvent;
public:
	ImmuneHistory(Host * hostPtr, int64_t const locusNumber, double infectionTimesToImmune);
	
	void gainImmunity(GenePtr genePtr);
	void gainAlleleImmunity(GenePtr genePtr, bool writeToDatabase,Database & db,zppdb::Table<AlleleImmunityRow> & table);
    void setAlleleLossEvent(int64_t locusId, int64_t AlleleId,double lossrate);
    double checkGeneImmunity(GenePtr genePtr);
	void gainGeneralImmunity();
    double checkGeneralImmunity(double a, double b);
    double checkGeneralImmunity(std::vector<double> params, double & immuneRate);
    void loseImmunity(GenePtr genePtr);
    void loseAlleleImmune(int64_t & locusId,int64_t & AlleleId);
	bool isImmune(GenePtr genePtr);
	
	void prepareToDie();
	
	void write(Database & db, Table<ImmunityRow> & table);
	void write(int64_t transmissionId, Database & db, Table<TransmissionImmunityRow> & table);
	
    // Vector of maps such that
    // immuneAlleles[i] is a map: locus number -> immunity level;
    // immuneAlleles[i][j] is the immunity level at locus i for allele j
    std::vector<std::unordered_map<int64_t,int64_t>> immuneAlleles;
    
	std::unordered_set<GenePtr> genes;
	std::unordered_map<GenePtr, std::unique_ptr<ImmunityLossEvent>> lossEvents;
    int64_t geneLocToHash(int64_t locusId, int64_t AlleleId);
	std::unordered_map<int64_t, std::unique_ptr<AlleleImmuneLossEvent>> alleleLossEvents;
    void updateAlleleLossRate(int64_t &  locusId, int64_t &  AlleleId,double newRate);
private:
	Host * hostPtr;
    int64_t const locusNumber;
    int64_t infectedTimes = 0;
    double infectionTimesToImmune;
};

#endif /* defined(__malariamodel__ImmuneHistory__) */
