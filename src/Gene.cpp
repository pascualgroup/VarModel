//
//  Gene.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 4/4/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "Gene.h"
#include "zppsim_random.hpp"
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;
using namespace zppsim;
using namespace zppdb;

Gene::Gene(
           int64_t id, double transmissibility, double immunityLossRate, int64_t const source,
           bool functionality, std::vector<int64_t> knownAlleles,bool writeToDatabaseGene, bool writeToDatabaseLoci, Database & db, Table<GeneRow> & GeneTable,Table<LociRow> & LociTable
) :
	id(id),
	transmissibility(transmissibility),
	immunityLossRate(immunityLossRate), source(source),functionality(functionality),
    Alleles(knownAlleles),recorded(writeToDatabaseGene)
{
	if(writeToDatabaseGene) {
		GeneRow row;
		row.geneId = id;
		row.transmissibility = transmissibility;
		row.immunityLossRate = immunityLossRate;
		row.source = source;
        row.functionality = functionality;
		db.insert(GeneTable, row);
	}
    if(writeToDatabaseLoci) {
        for(int64_t i = 0; i < Alleles.size(); i++) {
            LociRow row;
            row.geneId = id;
            row.alleleIndex = i;
            row.alleleId = Alleles[i];
            db.insert(LociTable,row);
        }
    }
}

void Gene::writeToDatabaseGene(Database & db, Table<GeneRow> & GeneTable,Table<LociRow> & LociTable) {
    if (! recorded){
        GeneRow row;
        row.geneId = id;
        row.transmissibility = transmissibility;
        row.immunityLossRate = immunityLossRate;
        row.source = source;
        row.functionality = functionality;
        db.insert(GeneTable, row);
        for(int64_t i = 0; i < Alleles.size(); i++) {
            LociRow row;
            row.geneId = id;
            row.alleleIndex = i;
            row.alleleId = Alleles[i];
            db.insert(LociTable,row);
        }
        recorded = true;
    }
    
}


std::string Gene::toString()
{
	stringstream ss;
	ss << id;
	return ss.str();
}
