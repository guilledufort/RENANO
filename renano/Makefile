OS_NAME := $(shell uname -s | tr A-Z a-z)

DEBUG_FLAGS = -pg
COMMON_FLAGS = -fopenmp -std=c++11 -lz -O3 -march=native -fstrict-aliasing -ffast-math -lz
OP_FLAGS = -fomit-frame-pointer
W_FLAGS = -Wall -Wextra

ifeq ($(OS_NAME),linux)
	CXX = g++
else
	CXX = g++-9
endif

all: renano

renano: 
		$(CXX) $(W_FLAGS) $(COMMON_FLAGS) $(OP_FLAGS) main.cpp enc_algs.cpp compressor.cpp alignments.cpp paf_process.cpp -o renano

debug: 
		$(CXX) $(W_FLAGS) $(COMMON_FLAGS) $(DEBUG_FLAGS) main.cpp enc_algs.cpp compressor.cpp alignments.cpp paf_process.cpp -o renano_debug

clean:
		rm -f renano *.o
