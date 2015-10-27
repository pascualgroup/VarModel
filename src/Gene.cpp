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
           int64_t id, double transmissibility, double immunityLossRate, double clinicalImmunityLossRate,
           bool functionality, std::vector<int64_t> knownAlleles,bool writeToDatabaseGene, bool writeToDatabaseLoci, Database & db, Table<GeneRow> & GeneTable,Table<LociRow> & LociTable
) :
	id(id),
	transmissibility(transmissibility),
	immunityLossRate(immunityLossRate), clinicalImmunityLossRate(clinicalImmunityLossRate),functionality(functionality),
    Alleles(knownAlleles)
{
	if(writeToDatabaseGene) {
		GeneRow row;
		row.geneId = id;
		row.transmissibility = transmissibility;
		row.immunityLossRate = immunityLossRate;
		row.clinicalImmunityLossRate = clinicalImmunityLossRate;
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

std::string Gene::toString()
{
	stringstream ss;
	ss << id;
	return ss.str();
}
