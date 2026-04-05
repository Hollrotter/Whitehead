#include "Structure.hpp"

int main()
{
    double L = 1.5; // Length of the plate in m
    double D = 0.5; // Height of the plate in m
    double d = 0.2; // Diameter of the hole in m
    double h = 0.15; // Distance between hole and border in m (redundant)
    double t = 0.01; // Thickness of the plate in m
    double delta = (D-d)/6; // Location of refined mesh

    double E = 2.1e11; // Youngs modulus
    double nu = 0.3; // Poissons ratio

    double P = 20000; // Load (pressure) in N/m

    size_t n1 = 10;
    size_t n2 = 10;
    size_t n3 = 15;
    size_t n4 = 10;
    size_t n5 = 20;

    Point p0(0, 0);
    Point p1(d/2, 0);
    Point p2(sqrt(2)/2*d/2, sqrt(2)/2*d/2);
    Point p3(0, d/2);
    Point p4(d/2+delta, 0);
    Point p5(sqrt(2)/2*(d/2+delta), sqrt(2)/2*(d/2+delta));
    Point p6(0, d/2+delta);
    Point p7(D/2, 0);
    Point p8(D/2, D/2);
    Point p9(0, D/2);
    Point p10(L/2, 0);
    Point p11(L/2, D/2);

    Lagrange::CurveInterpolant chi1 (p2,  p1,  p0, d/2, n1);
    Lagrange::CurveInterpolant chi2 (p3,  p2,  p0, d/2, n2);
    Lagrange::CurveInterpolant chi3 (p1,  p4,  n3);
    Lagrange::CurveInterpolant chi4 (p2,  p5,  n3);
    Lagrange::CurveInterpolant chi5 (p3,  p6,  n3);
    Lagrange::CurveInterpolant chi6 (p5,  p4,  p0, d/2+delta, n1);
    Lagrange::CurveInterpolant chi7 (p6,  p5,  p0, d/2+delta, n2);
    Lagrange::CurveInterpolant chi8 (p4,  p7,  n4);
    Lagrange::CurveInterpolant chi9 (p6,  p9,  n4);
    Lagrange::CurveInterpolant chi10(p8,  p7,  n1);
    Lagrange::CurveInterpolant chi11(p9,  p8,  n2);
    Lagrange::CurveInterpolant chi12(p7,  p10, n5);
    Lagrange::CurveInterpolant chi13(p11, p10, n1);
    Lagrange::CurveInterpolant chi14(p8,  p11, n5);
    Lagrange::CurveInterpolant chi15(p5,  p8,  n4);

    Membrane m1({&chi1,  &chi3,  &chi6,  &chi4 });
    Membrane m2({&chi2,  &chi4,  &chi7,  &chi5 });
    Membrane m3({&chi6,  &chi8,  &chi10, &chi15});
    Membrane m4({&chi7,  &chi15, &chi11, &chi9 });
    Membrane m5({&chi10, &chi12, &chi13, &chi14});

    m1.boundary(Field::n22, &chi1,  BC::Dirichlet);
    m1.boundary(Field::n12, &chi1,  BC::Dirichlet);
    m1.boundary(Field::v2,  &chi3,  BC::Neumann);
    m1.boundary(Field::v1,  &chi3,  BC::Dirichlet);

    m2.boundary(Field::n22, &chi2,  BC::Dirichlet);
    m2.boundary(Field::n12, &chi2,  BC::Dirichlet);
    m2.boundary(Field::v2,  &chi5,  BC::Neumann);
    m2.boundary(Field::v1,  &chi5,  BC::Dirichlet);

    m3.boundary(Field::v2,  &chi8,  BC::Neumann);
    m3.boundary(Field::v1,  &chi8,  BC::Dirichlet);
    m3.boundary(Field::z,   &chi8,  BC::Neumann);

    m4.boundary(Field::v2,  &chi9,  BC::Neumann);
    m4.boundary(Field::v1,  &chi9,  BC::Dirichlet);
    m4.boundary(Field::n22, &chi11, BC::Dirichlet);
    m4.boundary(Field::n12, &chi11, BC::Dirichlet);

    m5.boundary(Field::v2,  &chi12, BC::Neumann);
    m5.boundary(Field::v1,  &chi12, BC::Dirichlet);
    m5.boundary(Field::n22, &chi13, BC::Dirichlet, P);
    m5.boundary(Field::n12, &chi13, BC::Dirichlet);
    m5.boundary(Field::n12, &chi14, BC::Dirichlet);
    m5.boundary(Field::n11, &chi14, BC::Dirichlet);

    Structure s({&m1, &m2, &m3, &m4, &m5});
    s.setgamma0(1);
    s.youngsModulus(E*t);
    s.poissonsRatio(nu);
    s.setIterations(500);
    
    s.planeStrain();

    s.output("plot/Data/Kirsch_v1",  Field::v1);
    s.output("plot/Data/Kirsch_v2",  Field::v2);
    s.output("plot/Data/Kirsch_n11", Field::n11);
    s.output("plot/Data/Kirsch_n12", Field::n12);
    s.output("plot/Data/Kirsch_n22", Field::n22);
    s.principalStresses("plot/Data/Kirsch_sx", "plot/Data/Kirsch_s1");
    s.principalStrains("plot/Data/Kirsch_vx", "plot/Data/Kirsch_gx", "plot/Data/Kirsch_g1");
}