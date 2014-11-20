#!/usr/bin/env python

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'zppjson', 'scripts'))
import zppjson_doxygen
sys.path.append(os.path.join(os.path.dirname(__file__), 'zppdb', 'scripts'))
import zppdb_doxygen

if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	
	def jsonOutFilename(filename):
		return '{0}_jsontypes.h'.format(os.path.splitext(filename)[0])
	
	def dbOutFilename(filename):
		return '{0}_dbtypes.h'.format(os.path.splitext(filename)[0])
	
	jsonFilenames = [
		'src/SimParameters.h'
	]
	
	dbFilenames = [
		'src/DatabaseTypes.h'
	]
	
	for filename in jsonFilenames:
		outFile = open(jsonOutFilename(filename), 'w')
		zppjson_doxygen.processFile(filename, outFile)
		outFile.close()
	
	for filename in dbFilenames:
		outFile = open(dbOutFilename(filename), 'w')
		zppdb_doxygen.processFile(filename, outFile)
		outFile.close()
	
	os.system('doxygen')
	for filename in jsonFilenames:
		os.remove(jsonOutFilename(filename))
	
	for filename in dbFilenames:
		os.remove(dbOutFilename(filename))

