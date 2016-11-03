#!/usr/bin/env python

import os
import sys

EXEC_NAME = sys.argv[1]

# Compilers in preference order: tuple (cCompiler, cppCompiler, flags)
COMPILERS = [
#Intel compiler does not currently work
#	('icc', 'icpc', ''),
	('clang', 'clang++', '-stdlib=libc++'),
	('gcc', 'g++', '')
]

def execPresent(execName):
	return os.system('which {0}'.format(execName)) == 0

if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	
	for compiler in COMPILERS:
		if execPresent(compiler[0]) and execPresent(compiler[1]):
			cCompiler, cppCompiler, flags = compiler
			break
	
	includeDirs = [
		'zppjson/src',
		'zppdb/src',
		'zppsim/src',
		'preprocessor/include',
		'libjson',
		'/apps/amd64/libraries/sqlite/3.8.11.1/include'
	]
	
	srcDirs = [
		'src',
		'zppjson/src',
		'zppdb/src',
		'zppsim/src'
	]
	
	os.system('mkdir -p bin')
	os.system('{0} -O3 -c libjson/json.c -o bin/libjson.o'.format(cCompiler))
	os.system('{0} -O3 -std=c++11 {1} {2} bin/libjson.o {3} -o bin/{4} -lsqlite3'.format(
		cppCompiler,
		flags,
		' '.join(['-I{0}'.format(os.path.expanduser(x)) for x in includeDirs]),
		' '.join(['{0}/*.cpp'.format(x) for x in srcDirs]),
		EXEC_NAME
	))
	os.system('rm bin/libjson.o')
