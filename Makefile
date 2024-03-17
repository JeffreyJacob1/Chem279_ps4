

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++11 

factorial: factorial.h
	$(CPP) $(CPPFLAGS) -c factorial.cpp

molecule: molecule.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp 

CNDO: CNDO.h
	$(CPP) $(CPPFLAGS) -c CNDO.cpp

test: factorial molecule CNDO
	$(CPP) $(CPPFLAGS) -o test test.cpp factorial.o molecule.o CNDO.o -larmadillo



all: test 

clean:
	rm -f *.o test 

	