CC=CC
VECT_FLAGS= -fopenmp -ffast-math -march=knl -mavx512f -mavx512pf -mavx512er -mavx512cd -fopt-info-vec-optimized=vecreport.dat
COMMON_FLAGS=-g -O3 -std=c++0x -o bs-solctra-seq -Wall
SOURCE=solctra_multinode.h solctra_multinode.cpp main_multinode.cpp utils.h utils.cpp
RPT_FLAGS=
FP_FLAGS=

multinode: $(SOURCE)
	$(CC) -O3 -std=c++11 $(VECT_FLAGS) -o bs-solctra-multinode $(SOURCE)
	cp bs-solctra-multinode results
	rm bs-solctra-multinode;

clean:
	rm results/bs-solctra-multinode;
	rm vecreport.dat;
