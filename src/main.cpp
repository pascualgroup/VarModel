#include <iostream>
#include <sstream>

#include "Database.hpp"
#include "SimParameters.h"
#include "Simulation.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>

using namespace std;
using namespace zppdata;
using namespace boost::property_tree;

int main(int argc, char **argv) {
	if(argc < 2) {
		cerr << "Usage: " << argv[0] << " <params-filename>" << endl;
		return 1;
	}
	
	// Load parameters from JSON file
	ifstream paramsStream(argv[1]);
	ptree paramsPtree;
	read_json(paramsStream, paramsPtree);
	SimParameters params;
	params.loadFromPtree(paramsPtree);
	
	// Generate random seed if set to 0
	if(params.randomSeed == 0) {
		random_device rd;
		uniform_int_distribution<uint32_t> ud(0, std::numeric_limits<uint32_t>::max());
		params.randomSeed = ud(rd);
	}
	
	cerr << "Loaded parameters:" << endl;
	write_json(cerr, params.dumpToPtree());
	
	Database db(params.dbFilename);
	
	unique_ptr<Simulation> simPtr(new Simulation(&params, &db));
	simPtr->run();
}
