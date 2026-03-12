# -*- Makefile -*-

STRING = lib/String.o lib/StringLinear.o lib/StringNonlinear.o lib/StringStructuralMatrix.o
MEMBRANE = lib/Membrane.o lib/MembraneBoundary.o lib/MembraneDerivative.o lib/MembraneLinear.o \
lib/MembraneSemilinear.o lib/MembraneNonlinear.o lib/MembraneStructuralMatrix.o lib/PlaneStrain.o
STRUCTURE = lib/Structure.o lib/StructureLinear.o lib/StructureSemilinear.o lib/StructureNonlinear.o \
lib/StructurePlaneStrain.o
CHEBYSHEV = lib/Chebyshev.o lib/ChebyshevFFT.o lib/ChebyshevPolynomial.o
LAGRANGE = lib/Lagrange.o
METRIC = lib/Christoffel.o lib/Metric.o lib/MetricCo.o lib/Jacobian.o
SPLINE = lib/B_Spline.o lib/setQR.o lib/Splinefit.o lib/SplinefitDiff.o
DVM = lib/Camber.o lib/DVM.o lib/DVMAerodynamicMatrix.o
VLM = lib/VLM.o lib/VLMAerodynamicMatrix.o lib/Vortex.o
AIRFOIL = lib/Airfoil.o lib/AirfoilAerodynamicMatrix.o
WING = lib/Wing.o lib/WingAerodynamicMatrix.o lib/WingBoundary.o
MISC = lib/misc.o lib/fastgl.o
OBJS = $(STRING) $(MEMBRANE) $(STRUCTURE) $(CHEBYSHEV) $(LAGRANGE) $(METRIC) $(SPLINE) $(DVM) $(VLM) \
$(AIRFOIL) $(WING) $(MISC)
BINS = $(OBJS) lib/libWhitehead.a
OPTIONS = g++ -Ofast -c
INCLUDE = -I./include -fopenmp -std=c++23
LDLIBS = -llibarmadillo

all: $(BINS)

lib/Airfoil.o: src/Airfoil.cpp
	$(OPTIONS) src/Airfoil.cpp $(INCLUDE) $(LDLIBS) -o lib/Airfoil.o

lib/AirfoilAerodynamicMatrix.o: src/AirfoilAerodynamicMatrix.cpp
	$(OPTIONS) src/AirfoilAerodynamicMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/AirfoilAerodynamicMatrix.o

lib/B_Spline.o: src/B_Spline.cpp
	$(OPTIONS) src/B_Spline.cpp $(INCLUDE) $(LDLIBS) -o lib/B_Spline.o

lib/Camber.o: src/Camber.cpp
	$(OPTIONS) src/Camber.cpp $(INCLUDE) $(LDLIBS) -o lib/Camber.o

lib/Chebyshev.o: src/Chebyshev.cpp
	$(OPTIONS) src/Chebyshev.cpp $(INCLUDE) $(LDLIBS) -o lib/Chebyshev.o

lib/ChebyshevFFT.o: src/ChebyshevFFT.cpp
	$(OPTIONS) src/ChebyshevFFT.cpp $(INCLUDE) $(LDLIBS) -o lib/ChebyshevFFT.o

lib/ChebyshevPolynomial.o: src/ChebyshevPolynomial.cpp
	$(OPTIONS) src/ChebyshevPolynomial.cpp $(INCLUDE) $(LDLIBS) -o lib/ChebyshevPolynomial.o

lib/Christoffel.o: src/Christoffel.cpp
	$(OPTIONS) src/Christoffel.cpp $(INCLUDE) $(LDLIBS) -o lib/Christoffel.o

lib/DVM.o: src/DVM.cpp
	$(OPTIONS) src/DVM.cpp $(INCLUDE) $(LDLIBS) -o lib/DVM.o

lib/DVMAerodynamicMatrix.o: src/DVMAerodynamicMatrix.cpp
	$(OPTIONS) src/DVMAerodynamicMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/DVMAerodynamicMatrix.o

lib/fastgl.o: src/fastgl.cpp
	$(OPTIONS) src/fastgl.cpp $(INCLUDE) $(LDLIBS) -o lib/fastgl.o

lib/Jacobian.o: src/Jacobian.cpp
	$(OPTIONS) src/Jacobian.cpp $(INCLUDE) $(LDLIBS) -o lib/Jacobian.o

lib/Lagrange.o: src/Lagrange.cpp
	$(OPTIONS) src/Lagrange.cpp $(INCLUDE) $(LDLIBS) -o lib/Lagrange.o

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

lib/misc.o: src/misc.cpp
	$(OPTIONS) src/misc.cpp $(INCLUDE) $(LDLIBS) -o lib/misc.o

lib/PlaneStrain.o: src/PlaneStrain.cpp
	$(OPTIONS) src/PlaneStrain.cpp $(INCLUDE) $(LDLIBS) -o lib/PlaneStrain.o

lib/setQR.o: src/setQR.cpp
	$(OPTIONS) src/setQR.cpp $(INCLUDE) $(LDLIBS) -o lib/setQR.o

lib/Splinefit.o: src/Splinefit.cpp
	$(OPTIONS) src/Splinefit.cpp $(INCLUDE) $(LDLIBS) -o lib/Splinefit.o

lib/SplinefitDiff.o: src/SplinefitDiff.cpp
	$(OPTIONS) src/SplinefitDiff.cpp $(INCLUDE) $(LDLIBS) -o lib/SplinefitDiff.o

lib/String.o: src/String.cpp
	$(OPTIONS) src/String.cpp $(INCLUDE) $(LDLIBS) -o lib/String.o

lib/StringLinear.o: src/StringLinear.cpp
	$(OPTIONS) src/StringLinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StringLinear.o

lib/StringNonlinear.o: src/StringNonlinear.cpp
	$(OPTIONS) src/StringNonlinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StringNonlinear.o

lib/StringStructuralMatrix.o: src/StringStructuralMatrix.cpp
	$(OPTIONS) src/StringStructuralMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/StringStructuralMatrix.o

lib/Structure.o: src/Structure.cpp
	$(OPTIONS) src/Structure.cpp $(INCLUDE) $(LDLIBS) -o lib/Structure.o
	
lib/StructureLinear.o: src/StructureLinear.cpp
	$(OPTIONS) src/StructureLinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StructureLinear.o

lib/StructureSemilinear.o: src/StructureSemilinear.cpp
	$(OPTIONS) src/StructureSemilinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StructureSemilinear.o

lib/StructureNonlinear.o: src/StructureNonlinear.cpp
	$(OPTIONS) src/StructureNonlinear.cpp $(INCLUDE) $(LDLIBS) -o lib/StructureNonlinear.o

lib/StructurePlaneStrain.o: src/StructurePlaneStrain.cpp
	$(OPTIONS) src/StructurePlaneStrain.cpp $(INCLUDE) $(LDLIBS) -o lib/StructurePlaneStrain.o

lib/VLM.o: src/VLM.cpp
	$(OPTIONS) src/VLM.cpp $(INCLUDE) $(LDLIBS) -o lib/VLM.o

lib/VLMAerodynamicMatrix.o: src/VLMAerodynamicMatrix.cpp
	$(OPTIONS) src/VLMAerodynamicMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/VLMAerodynamicMatrix.o

lib/Vortex.o: src/Vortex.cpp
	$(OPTIONS) src/Vortex.cpp $(INCLUDE) $(LDLIBS) -o lib/Vortex.o

lib/Wing.o: src/Wing.cpp
	$(OPTIONS) src/Wing.cpp $(INCLUDE) $(LDLIBS) -o lib/Wing.o

lib/WingAerodynamicMatrix.o: src/WingAerodynamicMatrix.cpp
	$(OPTIONS) src/WingAerodynamicMatrix.cpp $(INCLUDE) $(LDLIBS) -o lib/WingAerodynamicMatrix.o

lib/WingBoundary.o: src/WingBoundary.cpp
	$(OPTIONS) src/WingBoundary.cpp $(INCLUDE) $(LDLIBS) -o lib/WingBoundary.o

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

cleanMisc:
	-rm $(MISC)