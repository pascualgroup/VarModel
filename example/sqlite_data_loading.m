% Get mksqlite here:
% http://sourceforge.net/projects/mksqlite/

% Create connection to SQLite database
dbid = mksqlite('open', 'test_database.sqlite')

% Retrieve all rows from a database table
results = mksqlite(dbid, 'SELECT * FROM hostSamples');

% Retrieve a subset of rows, including their unique rowids,
% from a database table using a conditional statement
results = mksqlite(dbid, 'SELECT rowid,* FROM hostSamples WHERE time > 0.5');

% Update time in row 2 (1-indexed like Matlab)
mksqlite(dbid, 'BEGIN TRANSACTION')
mksqlite(dbid, 'UPDATE hostSamples SET time = 0.4 WHERE rowid = 2')
mksqlite(dbid, 'COMMIT')

% Update time in row 3 using a variable
t = 0.3;
rowid = 3;
mksqlite(dbid, 'BEGIN TRANSACTION')
mksqlite(dbid, 'UPDATE hostSamples SET time = ? WHERE rowid = ?', t, rowid);
mksqlite(dbid, 'COMMIT')

% Close connection cleanly
mksqlite(dbid, 'close')
