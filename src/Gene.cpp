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
    bool writeToDatabase, Database & db, Table<GeneRow> & table
) :
    id(id),
    transmissibility(transmissibility),
    immunityLossRate(immunityLossRate), clinicalImmunityLossRate(clinicalImmunityLossRate)
{
    if(writeToDatabase) {
        GeneRow row;
        row.geneId = id;
        row.transmissibility = transmissibility;
        row.immunityLossRate = immunityLossRate;
        row.clinicalImmunityLossRate = clinicalImmunityLossRate;
        
        db.insert(table, row);
    }
}

std::string Gene::toString()
{
    stringstream ss;
    ss << id;
    return ss.str();
}
