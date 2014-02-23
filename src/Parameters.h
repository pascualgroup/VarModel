#ifndef __malariamodel__Parameters__
#define __malariamodel__Parameters__

#include <map>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include "Database.h"
#include "strutils.h"

#define LOAD_VALUE(x) { loadValue(x, #x); };
#define LOAD_MAP(x) { loadMap(x, #x); };
#define LOAD_VECTOR(x, size) { loadVector(x, #x, size); };
#define LOAD_MATRIX(x, nrows, ncols) { loadMatrix(x, #x, nrows, ncols); };

class Parameters
{
public:
	Parameters(std::istream & paramStream, Database & db, char const * const tableName):
		dbPtr(&db),
		table(tableName)
	{
		boost::property_tree::json_parser::read_json(paramStream, paramPTree);
		
		table.addColumn("param_name", DBType::TEXT);
		table.addColumn("param_value", DBType::TEXT);
		
		db.createTable(table);
	}

protected:
	Database * dbPtr;
	boost::property_tree::ptree paramPTree;
	DBTable table;
	
	template<typename T>
	void loadValue(T & var, char const * paramName)
	{
		var = paramPTree.get<T>(paramName);
		writeToDatabase(paramName, var);
	}
	
	template<typename T>
	void loadMap(std::map<std::string, T> var, char const * paramName)
	{
		auto tlItr = paramPTree.find(paramName);
		for(auto itr = tlItr->second.begin(); itr != tlItr->second.end(); itr++) {
			var[itr->first] = tlItr->second.get<T>(itr->first.c_str());
		}
		writeToDatabase(paramName, var);
	}
	
	template<typename T>
	void loadVector(std::vector<T> & vec, char const * paramName, uint32_t size)
	{
		vec.clear();
		try {
			// If specified as single double value, replicate to whole array
			T value = paramPTree.get<T>(paramName);
			for(uint32_t i = 0; i < size; i++) {
				vec.push_back(value);
			}
		}
		catch(boost::property_tree::ptree_error e) {
			// Otherwise load vector and make sure it's the right size
			auto tlItr = paramPTree.find(paramName);
			for(auto itr = tlItr->second.begin(); itr != tlItr->second.end(); itr++) {
				vec.push_back(boost::lexical_cast<double>(itr->second.data()));
			}
			assert(vec.size() == size);
		}
		writeToDatabase(paramName, vec);
	}
	
	template<typename T>
	void loadMatrix(std::vector<std::vector<T>> & mat, char const * paramName, uint32_t nrows, uint32_t ncols)
	{
		mat.clear();
		try {
			// If specified as a single double value, replicate to whole matrix
			T value = paramPTree.get<T>(paramName);
			mat = std::vector<std::vector<T>>(nrows, std::vector<T>(ncols, value));
		}
		catch(boost::property_tree::ptree_error e) {
			// Otherwise load matrix and make sure it's the right size
			auto tlItr = paramPTree.find(paramName);
			uint32_t i = 0;
			for(auto rowItr = tlItr->second.begin(); rowItr != tlItr->second.end(); rowItr++) {
				mat.push_back(std::vector<T>());
				for(auto colItr = rowItr->second.begin(); colItr != rowItr->second.end(); colItr++) {
					mat[i].push_back(boost::lexical_cast<double>(colItr->second.data()));
				}
				assert(mat[i].size() == ncols);
				i++;
			}
			assert(mat.size() == nrows);
		}
		writeToDatabase(paramName, mat);
	}
	
	template<typename T>
	void writeToDatabase(std::string paramName, T const & value)
	{
		DBRow row;
		row.set("param_name", paramName);
		row.set("param_value", toJson(value));
		dbPtr->insert(table, row);
	}
	
	std::string toJson(std::string const & value)
	{
		return strprintf("\"%s\"", escapeJsonString(value).c_str());
	}
	
	std::string toJson(bool const value)
	{
		return value ? "true" : "false";
	}
	
	std::string toJson(double const value)
	{
		return strprintf("%g", value);
	}
	
	std::string toJson(uint32_t const value)
	{
		return strprintf("%u", value);
	}
	
	template<typename T>
	std::string toJson(std::map<std::string, T> & value)
	{
		std::string jsonStr;
		jsonStr.append("{");
		std::vector<std::string> pieces;
		for(auto itr = value.begin(); itr != value.end(); itr++) {
			pieces.push_back(strprintf("\"%s\" : %s", escapeJsonString(itr->first).c_str(), toJson(itr->second)));
		}
		jsonStr.append(strjoin(pieces, ", "));
		jsonStr.append("}");
		
		return jsonStr;
	}
	
	template<typename T>
	std::string toJson(std::vector<T> const & value)
	{
		std::string jsonStr;
		jsonStr.append("[");
		std::vector<std::string> pieces;
		for(auto itr = value.begin(); itr != value.end(); itr++) {
			pieces.push_back(toJson(*itr));
		}
		jsonStr.append(strjoin(pieces, ", "));
		jsonStr.append("]");
		
		return jsonStr;
	}
	
	template<typename T>
	std::string toJson(std::vector<std::vector<T>> const & value)
	{
		std::vector<std::string> rowStrs;
		for(auto vecItr = value.begin(); vecItr != value.end(); vecItr++) {
			std::vector<std::string> dblStrs;
			for(auto itr = vecItr->begin(); itr != vecItr->end(); itr++) {
				dblStrs.push_back(toJson(*itr));
			}
			rowStrs.push_back(strprintf("[%s]", strjoin(dblStrs, ", ").c_str()));
		}
		std::string jsonStr = strprintf("[%s]", strjoin(rowStrs, ", ").c_str());
		
		return jsonStr;
	}
};

#endif
