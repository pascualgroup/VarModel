#include "Database.h"

#include <vector>
#include <iostream>
#include "strutils.h"

using namespace std;

DBTable::DBTable(std::string name) : name(name) {
}

std::string DBTable::getName() const {
	return name;
}

void DBTable::addColumn(std::string colName, DBType type) {
	size_t index = columns.size();
	columns.push_back(DBColumn(colName, type));
	indexMap[colName] = index;
}

bool DBTable::hasColumn(std::string colName) const {
	return indexMap.find(colName) != indexMap.end();
}

size_t DBTable::columnCount() const {
	return columns.size();
}

DBType DBTable::getColumnType(std::string colName) const {
	return getColumn(colName).type;
}

DBColumn DBTable::getColumn(std::string colName) const {
	return columns[indexMap.at(colName)];
}

DBColumn DBTable::getColumn(size_t index) const {
	return columns[index];
}

DBRow::DBRow() {
}

void DBRow::set_null(std::string const & key) {
	intMap.erase(key);
	dblMap.erase(key);
	strMap.erase(key);
}

void DBRow::set(std::string const & key, int64_t const & value) {
	intMap[key] = value;
}

bool DBRow::hasInteger(std::string const & key) const {
	return intMap.find(key) != intMap.end();
}

int64_t DBRow::getInteger(std::string const & key) const {
	return intMap.at(key);
}

void DBRow::set(std::string const & key, double const & value) {
	dblMap[key] = value;
}

bool DBRow::hasReal(std::string const & key) const {
	return dblMap.find(key) != dblMap.end();
}

double DBRow::getReal(std::string const & key) const {
	return dblMap.at(key);
}

void DBRow::set(std::string const & key, std::string const & value) {
	strMap[key] = value;
}

bool DBRow::hasText(std::string const & key) const {
	return strMap.find(key) != strMap.end();
}

std::string DBRow::getText(std::string const & key) const {
	return strMap.at(key);
}

bool DBRow::matches(DBTable const & table) const {
	for(auto & kv : intMap) {
		if(table.hasColumn(kv.first)) {
			if(table.getColumnType(kv.first) != DBType::INTEGER) {
				cerr << kv.first <<  "not integer" << endl;
				return false;
			}
		}
		else {
			cerr << kv.first << " not present" << endl;
			return false;
		}
	}
	
	for(auto & kv : dblMap) {
		if(table.hasColumn(kv.first)) {
			if(table.getColumnType(kv.first) != DBType::REAL) {
				cerr << kv.first << " not real" << endl;
				return false;
			}
		}
		else {
			cerr << kv.first << " not present" << endl;
			return false;
		}
	}
	
	for(auto & kv : strMap) {
		if(table.hasColumn(kv.first)) {
			if(table.getColumnType(kv.first) != DBType::TEXT) {
				cerr << kv.first << "not text" << endl;
				return false;
			}
		}
		else {
			cerr << kv.first << " not present" << endl;
			return false;
		}
	}
	
	return true;
}

Database::Database(bool enabled, string const & filename) {
	if(!enabled || filename == "") {
		db = NULL;
	}
	else {
		int result = sqlite3_open(filename.c_str(), &db);
		if(result) {
			throw std::runtime_error(string(sqlite3_errmsg(db)));
		}
	}
}

Database::~Database() {
	if(db != NULL) {
		sqlite3_close(db);
	}
}

void Database::beginTransaction() const {
	if(db == NULL) {
		return;
	}
	
	int result = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
	if(result) {
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
}

void Database::commit() const {
	if(db == NULL) {
		return;
	}
	
	int result = sqlite3_exec(db, "COMMIT", NULL, NULL, NULL);
	if(result) {
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
}
	
void Database::createTable(DBTable const & table) const {
	if(db == NULL) {
		return;
	}
	
	vector<string> columnStrs;
	for(size_t i = 0; i < table.columnCount(); i++) {
		DBColumn column = table.getColumn(i);
		switch(column.type) {
			case DBType::INTEGER:
				columnStrs.push_back(column.name + " INTEGER");
				break;
			case DBType::REAL:
				columnStrs.push_back(column.name + " REAL");
				break;
			case DBType::TEXT:
				columnStrs.push_back(column.name + " TEXT");
				break;
		}
	}
	string sqlCmd = strprintf(
		"CREATE TABLE %s (%s)",
		table.getName().c_str(),
		strjoin(columnStrs, ", ").c_str()
	);
	int result = sqlite3_exec(db, sqlCmd.c_str(), NULL, NULL, NULL);
	if(result) {
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
}

void Database::insert(DBTable const & table, DBRow const & row) const {
	if(!row.matches(table)) {
		throw std::runtime_error(
			strprintf("Row does not match specification for table %s.", table.getName().c_str())
		);
	}
	
	if(db == NULL) {
		return;
	}
	
	// Construct SQL statement with ? parameters
	string qMark("?");
	vector<string> qMarks;
	for(uint32_t i = 0; i < table.columnCount(); i++) {
		qMarks.push_back(qMark);
	}
	string sqlCmdStr = strprintf("INSERT INTO %s VALUES (%s)", table.getName().c_str(), strjoin(qMarks, ", ").c_str());
	char const * sqlCmd = sqlCmdStr.c_str();
	
	// Create prepared statement
	sqlite3_stmt *ppStmt;
	int result = sqlite3_prepare_v2(db, sqlCmd, -1, &ppStmt, NULL);
	if(result) {
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
	
	// Bind values to prepared statement parameters
	for(size_t i = 0; i < table.columnCount(); i++) {
		DBColumn col = table.getColumn(i);
		switch(col.type) {
			case DBType::INTEGER:
				if(row.hasInteger(col.name)) {
					result = sqlite3_bind_int64(ppStmt, int(i+1), row.getInteger(col.name));
				}
				else {
					result = sqlite3_bind_null(ppStmt, int(i+1));
				}
				break;
			case DBType::REAL:
				if(row.hasReal(col.name)) {
					result = sqlite3_bind_double(ppStmt, int(i+1), row.getReal(col.name));
				}
				else {
					result = sqlite3_bind_null(ppStmt, int(i+1));
				}
				break;
			case DBType::TEXT:
				if(row.hasText(col.name)) {
					string textStr = row.getText(col.name);
					char const * text = textStr.c_str();
					result = sqlite3_bind_text(ppStmt, int(i+1), text, -1, SQLITE_TRANSIENT);
				}
				else {
					result = sqlite3_bind_null(ppStmt, int(i+1));
				}
				break;
		}
		if(result) {
			throw std::runtime_error(string(sqlite3_errmsg(db)));
		}
	}
	
	// Execute prepared statement
	result = sqlite3_step(ppStmt);
	if(result != SQLITE_DONE) {
		sqlite3_finalize(ppStmt);
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
	
	sqlite3_finalize(ppStmt);
}

void Database::insert(DBTable const & table, std::vector<double> & row) const
{
	if(table.columnCount() != row.size()) {
		throw std::runtime_error(
			strprintf("Row does not match column count in table %s.", table.getName().c_str())
		);
	}
	
	if(db == NULL) {
		return;
	}
	
	// Construct SQL statement with ? parameters
	string qMark("?");
	vector<string> qMarks;
	for(uint32_t i = 0; i < table.columnCount(); i++) {
		qMarks.push_back(qMark);
	}
	string sqlCmdStr = strprintf("INSERT INTO %s VALUES (%s)", table.getName().c_str(), strjoin(qMarks, ", ").c_str());
	char const * sqlCmd = sqlCmdStr.c_str();
	
	// Create prepared statement
	sqlite3_stmt *ppStmt;
	int result = sqlite3_prepare_v2(db, sqlCmd, -1, &ppStmt, NULL);
	if(result) {
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
	
	// Bind values to prepared statement parameters
	for(size_t i = 0; i < table.columnCount(); i++) {
		result = sqlite3_bind_double(ppStmt, int(i+1), row[i]);
		if(result) {
			throw std::runtime_error(string(sqlite3_errmsg(db)));
		}
	}
	
	// Execute prepared statement
	result = sqlite3_step(ppStmt);
	if(result != SQLITE_DONE) {
		sqlite3_finalize(ppStmt);
		throw std::runtime_error(string(sqlite3_errmsg(db)));
	}
	
	sqlite3_finalize(ppStmt);
}
