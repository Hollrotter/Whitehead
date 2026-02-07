# -*- Makefile -*-

STRING = lib/String.o lib/StringLinear.o lib/StringNonlinear.o lib/StringStructuralMatrix.o
CHEBYSHEV = lib/Chebyshev.o
OBJS = $(STRING) $(CHEBYSHEV)
BINS = $(OBJS) lib/libWhitehead.a
OPTIONS = g++ -Ofast -c
INCLUDE = -I./include -fopenmp -std=c++23
LDLIBS = -llibarmadillo

all: $(BINS)

lib/Chebyshev.o: src/Chebyshev.cpp
	$(OPTIONS) src/Chebyshev.cpp $(INCLUDE) $(LDLIBS) -o lib/Chebyshev.o

lib/String.o: src/String.cpp
	$(OPTIONS) src/String.cpp $(INCLUDE) $(LDLIBS) -o lib/String.o

lib/StringLinear.o: src/StringLinear.cpp
	$(OPTIONS) src/StringLinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StringLinear.o

lib/StringNonlinear.o: src/StringNonlinear.cpp
	$(OPTIONS) src/StringNonlinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StringNonlinear.o

lib/StringStructuralMatrix.o: src/StringStructuralMatrix.cpp
	$(OPTIONS) src/StringStructuralMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/StringStructuralMatrix.o

lib/libWhitehead.a: $(OBJS)
	ar rcs lib/libWhitehead.a $(OBJS)

clean:
	-rm $(BINS)

cleanString:
	-rm $(STRING)

cleanChebyshev:
	-rm $(CHEBYSHEV)