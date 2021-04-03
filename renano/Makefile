OS_NAME := $(shell uname -s | tr A-Z a-z)

DEBUG_FLAGS = -pg
CFLAGS = -fopenmp -std=c++11 -O3 -march=native -fstrict-aliasing -ffast-math -fomit-frame-pointer -Wall
LDFLAGS  += -lz
INCLUDES=

ifeq ($(OS_NAME),linux)
	CXX = g++
else
	CXX = g++-9
endif

all: renano

renano: 
		$(CXX) $(CFLAGS) $(LDFLAGS) $(INCLUDES) main.cpp enc_algs.cpp compressor.cpp alignments.cpp paf_process.cpp -o renano

debug: 
		$(CXX) $(CFLAGS) $(LDFLAGS) $(INCLUDES) main.cpp enc_algs.cpp compressor.cpp alignments.cpp paf_process.cpp -o renano_debug

clean:
		rm -f renano *.o
