/**
	\mainpage malariamodel
	
	Code documentation for malaria var-gene evolution model.
	
	Useful sections:
	- Parameters: SimParameters_jsontypes.h
	- Database tables: DatabaseTypes_dbtypes.h: database table types
 */

#include <iostream>
#include <sstream>

#include "zppjson.hpp"
#include "zppdb.hpp"
#include "DatabaseTypes.h"
#include "SimParameters.h"
#include "Simulation.h"

#include <iostream>
#include <fstream>
#include <random>

using namespace std;
using namespace zppjson;
using namespace zppdb;

/**
	\brief Main entry point to program.
	
	Loads parameters, initializes random seed, initializes and runs simulation.
 */
int main(int argc, char **argv)
{
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " <params-filename>" << endl;
        return 1;
    }
    
    std::string paramsFilename(argv[1]);
    ifstream paramsStream(paramsFilename.c_str());
    if(!paramsStream) {
        cerr << "Parameters file " << paramsFilename << " cannot be opened.";
        return 1;
    }
    
    try {
        // Load parameters from JSON file
        JsonElement paramsJson = parseJson(paramsStream);
        paramsStream.close();
        
        SimParameters params;
        params.loadJson(paramsJson);
        
        // Generate random seed if set to 0
        if(!params.randomSeed.present() || (params.randomSeed == 0)) {
            random_device rd;
            uniform_int_distribution<uint32_t> ud(0, std::numeric_limits<uint32_t>::max());
            params.randomSeed = ud(rd);
        }
        
        stringstream paramsStrStream;
        cerr << "Loaded parameters:" << endl;
        params.printJsonToStream(cerr);
        params.printJsonToStream(paramsStrStream);
        
        // Set up database
        if(params.overwriteDatabase.present() && params.overwriteDatabase) {
            string dbFilename = params.dbFilename;
            remove(dbFilename.c_str());
            string journalFilename = dbFilename + "-journal";
            remove(journalFilename.c_str());
        }
        
        Database db(params.dbFilename);
        
        // Write parameters to meta table
        db.beginTransaction();
        Table<MetaRow> metaTable("meta");
        db.createTable(metaTable);
        MetaRow row;
        row.key = "parameters";
        row.value = paramsStrStream.str();
        db.insert(metaTable, row);
        db.commit();
        
        // Set up and run simulation
        unique_ptr<Simulation> simPtr(new Simulation(&params, &db));
        simPtr->run();
    }
    catch(runtime_error & e) {
        cerr << "An error occurred:" << endl;
        cerr << e.what() << endl;
        cerr << "Terminating execution." << endl;
        return 1;
    }
}
