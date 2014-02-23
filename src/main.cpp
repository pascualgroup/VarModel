#include <iostream>
#include <sstream>

#include "Database.h"
#include "SimParameters.h"
#include "Simulation.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>

using namespace std;
using namespace boost::property_tree;

int main(int argc, char **argv) {
	if(argc < 2) {
		cerr << "Usage: " << argv[0] << " <params-filename>" << endl;
		return 1;
	}
	
	Database db(true, "db.sqlite");
	
	ifstream paramsStream(argv[1]);
	SimParameters params(paramsStream, db);
	
	// Generate random seed if set to 0
	if(params.randomSeed == 0) {
		random_device rd;
		uniform_int_distribution<uint32_t> ud(0, std::numeric_limits<uint32_t>::max());
		params.randomSeed = ud(rd);
	}
	
	unique_ptr<Simulation> simulation(new Simulation(params, db));
	simulation->run();
}
