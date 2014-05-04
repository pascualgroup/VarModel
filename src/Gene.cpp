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
#include "zppdata_util.hpp"

using namespace std;
using namespace zppsim;
using namespace zppdata;

Gene::Gene(size_t id, double transmissibility, double immunityLossRate) :
	id(id), transmissibility(transmissibility), immunityLossRate(immunityLossRate)
{
}

std::string Gene::toString()
{
	return strprintf("g%u", id);
}
