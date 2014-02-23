all: build

build:
	c++ -O3 -std=c++11 -stdlib=libc++ -I. -lsqlite3 src/*.cpp -o malariamodel

clean:
	rm malariamodel
