//
//  Loci.h
//  malariamodel
//
//  Created by Qixin on 8/24/15.
//  Copyright (c) 2015 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__Loci__
#define __malariamodel__Loci__

#include <memory>
#include <vector>
#include <unordered_map>
#include "zppdb.hpp"
#include "zppjson.hpp"
#include "zppsim_random.hpp"
#include "DatabaseTypes.h"

//using namespace zppjson;

class Loci;
typedef std::shared_ptr<Loci> LociPtr;
typedef std::weak_ptr<Loci> LociPtrW;


class Loci
{
friend class Gene;
friend class Simulation;
public:
    int64_t const id;
    bool const functionality;
    //construct loci using new randomly generated alleles, with different allele numbers per locus
    Loci(zppsim::rng_t & rng, int64_t id,int64_t locusNumber, std::vector<int64_t> & alleleNumber, bool writeToDatabase, Database & db, Table<LociRow> & table);
    //construct loci using known alleles
    Loci(int64_t id,std::vector<int64_t> Alleles, bool functionality, bool writeToDatabase, Database & db, Table<LociRow> & table);
private:
    std::vector<int64_t> Alleles;
};


#endif /* defined(__malariamodel__Loci__) */
