# -*- Makefile -*-

CHEBYSHEV = lib/Chebyshev.o
OBJS = $(CHEBYSHEV)
BINS = $(OBJS) lib/libWhitehead.a
OPTIONS = g++ -Ofast -c
INCLUDE = -I./include -fopenmp -std=c++23
LDLIBS = -llibarmadillo

all: $(BINS)

lib/Chebyshev.o: src/Chebyshev.cpp
	$(OPTIONS) src/Chebyshev.cpp $(INCLUDE) $(LDLIBS) -o lib/Chebyshev.o

lib/libWhitehead.a: $(OBJS)
	ar rcs lib/libWhitehead.a $(OBJS)

clean:
	-rm $(BINS)

cleanChebyshev:
	-rm $(CHEBYSHEV)