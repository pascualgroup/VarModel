#ifndef __malariamodel__Database__
#define __malariamodel__Database__

#include <string>
#include <vector>
#include <unordered_map>
#include <sqlite3.h>

enum class DBType {
	INTEGER,
	REAL,
	TEXT
};

struct DBColumn {
	DBColumn(std::string name, DBType type) : name(name), type(type) { }
	std::string name;
	DBType type;
};

class DBTable {
public:
	DBTable(std::string name);
	std::string getName() const;
	void addColumn(std::string colName, DBType type);
	bool hasColumn(std::string colName) const;
	size_t columnCount() const;
	DBType getColumnType(std::string colName) const;
	DBColumn getColumn(std::string colName) const;
	DBColumn getColumn(size_t index) const;
private:
	std::string name;
	std::vector<DBColumn> columns;
	std::unordered_map<std::string, size_t> indexMap;
};

class DBRow {
public:
	DBRow();
	
	void set_null(std::string const & key);
	
	void set(std::string const & key, int64_t const & value);
	bool hasInteger(std::string const & key) const;
	int64_t getInteger(std::string const & key) const;
	
	void set(std::string const & key, double const & value);
	bool hasReal(std::string const & key) const;
	double getReal(std::string const & key) const;
	
	void set(std::string const & key, std::string const & value);
	bool hasText(std::string const & key) const;
	std::string getText(std::string const & key) const;
	
	bool matches(DBTable const & table) const;
private:
	std::unordered_map<std::string, int64_t> intMap;
	std::unordered_map<std::string, double> dblMap;
	std::unordered_map<std::string, std::string> strMap;
};

class Database {
public:
	Database(bool enabled, std::string const & filename);
	~Database();
	
	void beginTransaction() const;
	void commit() const;
	
	void createTable(DBTable const & table) const;
	void insert(DBTable const & table, DBRow const & row) const;
	void insert(DBTable const & table, std::vector<double> & row) const;
private:
	sqlite3 * db;
};

#endif
