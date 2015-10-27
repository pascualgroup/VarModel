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
    bool const functionality;
    std::vector<int64_t> Alleles;
	
	Gene(
		int64_t id, double transmissibility, double immunityLossRate, double clinicalImmunityLossRate,
		    bool functionality, std::vector<int64_t> knownAlleles,bool writeToDatabaseGene, bool writeToDatabaseLoci, Database & db, Table<GeneRow> & GeneTable,Table<LociRow> & LociTable
	);
	std::string toString();
};

#endif /* defined(__malariamodel__Gene__) */
