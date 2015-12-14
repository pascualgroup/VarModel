//
//  Gene.h
//  malariamodel
//
//  Created by Ed Baskerville on 4/4/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__Gene__
#define __malariamodel__Gene__

#include <memory>
#include <vector>
#include <unordered_map>
#include "zppdb.hpp"
#include "zppsim_random.hpp"
#include "DatabaseTypes.h"

class Gene;
typedef std::shared_ptr<Gene> GenePtr;
typedef std::weak_ptr<Gene> GenePtrW;

class Gene
{
public:
    int64_t const id;
    double const transmissibility;
    double const immunityLossRate;
    double const clinicalImmunityLossRate;
    
    Gene(
        int64_t id, double transmissibility, double immunityLossRate, double clinicalImmunityLossRate,
        bool writeToDatabase, Database & db, Table<GeneRow> & table
    );
    
    std::string toString();
};

#endif /* defined(__malariamodel__Gene__) */
