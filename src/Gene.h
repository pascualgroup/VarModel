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
#include "zppsim_random.hpp"

class Gene;
typedef std::shared_ptr<Gene> GenePtr;
typedef std::weak_ptr<Gene> GenePtrW;

class Gene
{
public:
	size_t const id;
	double const transmissibility;
	double const immunityLossRate;
	
	Gene(size_t id, double transmissibility, double immunityLossRate);
	
	std::string toString();
};

#endif /* defined(__malariamodel__Gene__) */
