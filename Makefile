ifneq ("$(wildcard /usr/bin/clang++)","")
	CPP = clang++
	EXTRA_FLAGS = -stdlib=libc++
else
	CPP = c++
	EXTRA_FLAGS = 
endif

all: build

build:
	$(CPP) -O3 -std=c++11 $(EXTRA_FLAGS) -lsqlite3 -Izppdata/src -Izppsim/src  src/*.cpp zppdata/src/*.cpp zppsim/src/*.cpp -o malariamodel

clean:
	rm malariamodel
