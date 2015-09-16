//
//  Loci.cpp
//  malariamodel
//
//  Created by Qixin on 8/24/15.
//  Copyright (c) 2015 Ed Baskerville. All rights reserved.
//

#include "Loci.h"
#include "zppsim_random.hpp"
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;
using namespace zppsim;
using namespace zppdb;

//construct loci using new randomly generated alleles
Loci::Loci(rng_t & rng, int64_t id,int64_t locusNumber, std::vector<int64_t> & alleleNumber,bool writeToDatabase, Database & db, Table<LociRow> & table) :
id(id),
functionality(true)
{
    Alleles.reserve(locusNumber);
    if(alleleNumber.size()==1) {
        std::uniform_int_distribution<int64_t> unif(0, alleleNumber[0]-1);
        for(int64_t i = 0; i < locusNumber; i++) {
            Alleles[i] = unif(rng);
        }
    }else{
        assert(alleleNumber.size()==locusNumber);
        for(int64_t i = 0; i < locusNumber; i++) {
            std::uniform_int_distribution<int64_t> unif(0, alleleNumber[i]-1);
            Alleles[i] = unif(rng);
        }
    }
    if(writeToDatabase) {
        for(int64_t i = 0; i < Alleles.size(); i++) {
            LociRow row;
            row.geneId = id;
            row.functionality = functionality;
            row.alleleIndex = i;
            row.alleleId = Alleles[i];
            db.insert(table,row);
        }
    }
}


//construct loci using known alleles
Loci::Loci(int64_t id,std::vector<int64_t> Alleles, bool functionality,bool writeToDatabase, Database & db, Table<LociRow> & table):
id(id),
Alleles(Alleles),
functionality(functionality)
{
    if(writeToDatabase) {
        for(int64_t i = 0; i < Alleles.size(); i++) {
            LociRow row;
            row.geneId = id;
            row.functionality = functionality;
            row.alleleIndex = i;
            row.alleleId = Alleles[i];
            db.insert(table,row);
        }
    }
}

