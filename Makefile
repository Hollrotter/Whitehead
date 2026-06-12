# -*- Makefile -*-

STRING = lib/String.o lib/StringLinear.o lib/StringNonlinear.o lib/StringStructuralMatrix.o
MEMBRANE = lib/Membrane.o lib/MembraneBoundary.o lib/MembraneDerivative.o lib/MembraneLinear.o \
lib/MembraneSemilinear.o lib/MembraneNonlinear.o lib/MembraneStructuralMatrix.o lib/PlaneStrain.o
STRUCTURE = lib/Structure.o lib/StructureLinear.o lib/StructureSemilinear.o lib/StructureNonlinear.o \
lib/StructurePlaneStrain.o
CHEBYSHEV = lib/Chebyshev.o lib/ChebyshevFFT.o lib/ChebyshevPolynomial.o
LAGRANGE = lib/Lagrange.o lib/CurveInterpolant.o
METRIC = lib/Christoffel.o lib/Metric.o lib/MetricCo.o lib/Jacobian.o
SPLINE = lib/B_Spline.o lib/setQR.o lib/Splinefit.o lib/SplinefitDiff.o
DVM = lib/Camber.o lib/DVM.o lib/DVMAerodynamicMatrix.o
VLM = lib/VLM.o lib/VLMAerodynamicMatrix.o lib/Vortex.o
AIRFOIL = lib/Airfoil.o lib/AirfoilAerodynamicMatrix.o lib/AirfoilKernel.o
WING = lib/Wing.o lib/WingAerodynamicMatrix.o lib/WingBoundary.o
AERODYNAMICS = lib/Aerodynamics.o lib/AerodynamicsBoundary.o lib/AerodynamicsLinear.o lib/AerodynamicsNonlinear.o
MISC = lib/misc.o lib/fastgl.o
OBJS = $(STRING) $(MEMBRANE) $(STRUCTURE) $(CHEBYSHEV) $(LAGRANGE) $(METRIC) $(SPLINE) $(DVM) $(VLM) \
$(AIRFOIL) $(WING) $(AERODYNAMICS) $(MISC)
BINS = $(OBJS) lib/libWhitehead.a
OPTIONS = g++ -Ofast -Wall -c
INCLUDE = -I./include -fopenmp -std=c++23

all: $(BINS)

lib/%.o: src/%.cpp
	$(OPTIONS) $(INCLUDE) -o $@ $^

lib/libWhitehead.a: $(OBJS)
	ar rcs lib/libWhitehead.a $(OBJS)

clean:
	-rm $(BINS)

cleanString:
	-rm $(STRING)

cleanMembrane:
	-rm $(MEMBRANE)

cleanStructure:
	-rm $(STRUCTURE)

cleanChebyshev:
	-rm $(CHEBYSHEV)

cleanMetric:
	-rm $(METRIC)

cleanLagrange:
	-rm $(LAGRANGE)

cleanSpline:
	-rm $(SPLINE)

cleanDVM:
	-rm $(DVM)

cleanVLM:
	-rm $(VLM)

cleanAirfoil:
	-rm $(AIRFOIL)

cleanWing:
	-rm $(WING)

cleanAerodynamics:
	-rm $(AERODYNAMICS)

cleanMisc:
	-rm $(MISC)