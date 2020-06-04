CXX = g++
CXXFLAGS = -g -std=c++11 -ffast-math -O2 -Wall -Wextra -pedantic-errors -Wold-style-cast

PROGRAMS = main tests

OBJECTS = deposition.o particle.o vector_operations.o inputreader.o fields_and_forces.o random.o

.DEFAULT_GOAL := default

.PHONY = all clean tests main docs

default: main

all: $(PROGRAMS)

main: main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS)  -o ../bin/AeroDep main.o $(OBJECTS) 
tests: ../tests/tests.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o ../tests/tests ../tests/tests.o $(OBJECTS) 

docs:
	doxygen ../docs/doxygen_config

clean: 
	rm -f *.o ../tests/*.o ../bin/AeroDep ../tests/tests
