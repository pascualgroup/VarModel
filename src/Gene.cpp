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

Gene::Gene(size_t id, double transmissibility, double immunityLossRate, double clinicalImmunityLossRate, DBTable * table) :
	id(id),
	transmissibility(transmissibility),
	immunityLossRate(immunityLossRate), clinicalImmunityLossRate(clinicalImmunityLossRate)
{
	if(table != nullptr) {
		DBRow row;
		row.setInteger("geneId", int64_t(id));
		row.setReal("transmissibility", transmissibility);
		row.setReal("immunityLossRate", immunityLossRate);
		row.setReal("clinicalImmunityLossRate", clinicalImmunityLossRate);
		table->insert(row);
	}
}

std::string Gene::toString()
{
	return strprintf("g%u", id);
}
