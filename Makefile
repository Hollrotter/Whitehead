# -*- Makefile -*-

STRING = lib/String.o lib/StringLinear.o lib/StringNonlinear.o lib/StringStructuralMatrix.o
MEMBRANE = lib/Membrane.o lib/MembraneBoundary.o lib/MembraneDerivative.o lib/MembraneLinear.o \
lib/MembraneSemilinear.o lib/MembraneNonlinear.o lib/MembraneStructuralMatrix.o lib/PlaneStrain.o
CHEBYSHEV = lib/Chebyshev.o
METRIC = lib/Christoffel.o lib/Metric.o lib/MetricCo.o lib/Jacobian.o
OBJS = $(STRING) $(MEMBRANE) $(CHEBYSHEV) $(METRIC)
BINS = $(OBJS) lib/libWhitehead.a
OPTIONS = g++ -Ofast -c
INCLUDE = -I./include -fopenmp -std=c++23
LDLIBS = -llibarmadillo

all: $(BINS)

lib/Chebyshev.o: src/Chebyshev.cpp
	$(OPTIONS) src/Chebyshev.cpp $(INCLUDE) $(LDLIBS) -o lib/Chebyshev.o

lib/Christoffel.o: src/Christoffel.cpp
	$(OPTIONS) src/Christoffel.cpp $(INCLUDE) $(LDLIBS) -o lib/Christoffel.o

lib/Jacobian.o: src/Jacobian.cpp
	$(OPTIONS) src/Jacobian.cpp $(INCLUDE) $(LDLIBS) -o lib/Jacobian.o

lib/Membrane.o: src/Membrane.cpp
	$(OPTIONS) src/Membrane.cpp $(INCLUDE) $(LDLIBS) -o lib/Membrane.o

lib/MembraneBoundary.o: src/MembraneBoundary.cpp
	$(OPTIONS) src/MembraneBoundary.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneBoundary.o

lib/MembraneDerivative.o: src/MembraneDerivative.cpp
	$(OPTIONS) src/MembraneDerivative.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneDerivative.o

lib/MembraneLinear.o: src/MembraneLinear.cpp
	$(OPTIONS) src/MembraneLinear.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneLinear.o

lib/MembraneSemilinear.o: src/MembraneSemilinear.cpp
	$(OPTIONS) src/MembraneSemilinear.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneSemilinear.o

lib/MembraneNonlinear.o: src/MembraneNonlinear.cpp
	$(OPTIONS) src/MembraneNonlinear.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneNonlinear.o

lib/MembraneStructuralMatrix.o: src/MembraneStructuralMatrix.cpp
	$(OPTIONS) src/MembraneStructuralMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/MembraneStructuralMatrix.o

lib/Metric.o: src/Metric.cpp
	$(OPTIONS) src/Metric.cpp $(INCLUDE) $(LDLIBS) -o lib/Metric.o

lib/MetricCo.o: src/MetricCo.cpp
	$(OPTIONS) src/MetricCo.cpp $(INCLUDE) $(LDLIBS) -o lib/MetricCo.o

lib/PlaneStrain.o: src/PlaneStrain.cpp
	$(OPTIONS) src/PlaneStrain.cpp $(INCLUDE) $(LDLIBS) -o lib/PlaneStrain.o

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