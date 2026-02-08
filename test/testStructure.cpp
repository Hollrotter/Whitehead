#include "Structure.hpp"

int main()
{
    switch(6)
    {
        case 0: // Rechteck
        {
            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;

            Point p1(0, 0);
            Point p2(1, 0);
            Point p3(1, 1);
            Point p4(0, 1);
            Point p5(2, 0);
            Point p6(2, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, n1);
            Lagrange::CurveInterpolant chi2(p2, p3, n2);
            Lagrange::CurveInterpolant chi3(p3, p4, n1);
            Lagrange::CurveInterpolant chi4(p1, p4, n2);

            Lagrange::CurveInterpolant chi5(p2, p5, n3);
            Lagrange::CurveInterpolant chi6(p5, p6, n2);
            Lagrange::CurveInterpolant chi7(p3, p6, n3);

            Membrane m1({&chi1, &chi2, &chi3, &chi4});

            m1.boundary(Field::z,   &chi1, BC::Neumann);
            m1.boundary(Field::v1,  &chi1, BC::Neumann);
            m1.boundary(Field::v2,  &chi1, BC::Dirichlet);

            m1.boundary(Field::z,   &chi3, BC::Dirichlet);
            m1.boundary(Field::n12, &chi3, BC::Dirichlet);
            m1.boundary(Field::n22, &chi3, BC::Dirichlet, 1);

            m1.boundary(Field::z,   &chi4, BC::Neumann);
            m1.boundary(Field::v1,  &chi4, BC::Dirichlet);
            m1.boundary(Field::v2,  &chi4, BC::Neumann);

            Membrane m2({&chi5, &chi6, &chi7, &chi2});

            m2.boundary(Field::z,   &chi5, BC::Neumann);
            m2.boundary(Field::v1,  &chi5, BC::Neumann);
            m2.boundary(Field::v2,  &chi5, BC::Dirichlet);

            m2.boundary(Field::z,   &chi6, BC::Dirichlet);
            m2.boundary(Field::n11, &chi6, BC::Dirichlet, 1);
            m2.boundary(Field::n12, &chi6, BC::Dirichlet);

            m2.boundary(Field::z,   &chi7, BC::Dirichlet);
            m2.boundary(Field::n12, &chi7, BC::Dirichlet);
            m2.boundary(Field::n22, &chi7, BC::Dirichlet, 1);

            Structure s({&m1, &m2});
            s.setgamma0(2.25);
            s.load(0.1);

            // s.inPlane11(10);
            // s.inPlane12(0);
            // s.inPlane22(10);

            s.setIterations(300);
            s.planeStrain();
            std::cout << "\n\n\n\n";
            s.linear();
            std::cout << "\n\n\n\n";
            s.semilinear();

            // m1.boundary(Field::n12, &chi3, BC::None);
            // m1.boundary(Field::n22, &chi3, BC::None);
            // m2.boundary(Field::n11, &chi6, BC::None);
            // m2.boundary(Field::n12, &chi6, BC::None);
            // m2.boundary(Field::n12, &chi7, BC::None);
            // m2.boundary(Field::n22, &chi7, BC::None);

            // m1.boundary(Field::v1, &chi3, BC::Dirichlet);
            // m1.boundary(Field::v2, &chi3, BC::Dirichlet);
            // m2.boundary(Field::v1, &chi6, BC::Dirichlet);
            // m2.boundary(Field::v2, &chi6, BC::Dirichlet);
            // m2.boundary(Field::v1, &chi7, BC::Dirichlet);
            // m2.boundary(Field::v2, &chi7, BC::Dirichlet);

            // s.nonlinear();

            s.output("plot/Data/Structure/v1", Field::v1);
            s.output("plot/Data/Structure/v2", Field::v2);
            s.output("plot/Data/Structure/z",  Field::z);
            s.principalStresses("plot/Data/Structure/kartesianStresses", "plot/Data/Structure/principalStresses");
            s.principalStrains("plot/Data/Structure/kartesianDeformations", "plot/Data/Structure/kartesianStrains", "plot/Data/Structure/principalStrains");
            break;
        }
        case 1: // Ring
        {
            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;

            Point p0(0, 0);
            Point p1(0, 1);
            Point p2(1, 0);
            Point p3(0, 2);
            Point p4(2, 0);
            Point p5(0, 3);
            Point p6(3, 0);

            Lagrange::CurveInterpolant chi1(p1, p3, n1);
            Lagrange::CurveInterpolant chi2(p1, p2, p0, 1, n2);
            Lagrange::CurveInterpolant chi3(p2, p4, n1);
            Lagrange::CurveInterpolant chi4(p3, p4, p0, 2, n2);
            Lagrange::CurveInterpolant chi5(p3, p5, n3);
            Lagrange::CurveInterpolant chi6(p5, p6, p0, 3, n2);
            Lagrange::CurveInterpolant chi7(p4, p6, n3);

            Membrane m1({&chi1, &chi2, &chi3, &chi4});
            Membrane m2({&chi5, &chi4, &chi7, &chi6});

            m1.boundary(Field::z,   &chi2, BC::Neumann);
            m1.boundary(Field::n11, &chi2, BC::Dirichlet);
            m1.boundary(Field::n12, &chi2, BC::Dirichlet);
            m1.boundary(Field::z,   &chi1, BC::Neumann);
            m1.boundary(Field::v2,  &chi1, BC::Dirichlet);
            m1.boundary(Field::v1,  &chi1, BC::Neumann);
            m1.boundary(Field::z,   &chi3, BC::Neumann);
            m1.boundary(Field::v2,  &chi3, BC::Dirichlet);
            m1.boundary(Field::v1,  &chi3, BC::Neumann);

            m2.boundary(Field::z,   &chi6, BC::Dirichlet);
            m2.boundary(Field::n12, &chi6, BC::Dirichlet);
            m2.boundary(Field::n11, &chi6, BC::Dirichlet, 1);
            m2.boundary(Field::z,   &chi5, BC::Neumann);
            m2.boundary(Field::v2,  &chi5, BC::Dirichlet);
            m2.boundary(Field::v1,  &chi5, BC::Neumann);
            m2.boundary(Field::z,   &chi7, BC::Neumann);
            m2.boundary(Field::v2,  &chi7, BC::Dirichlet);
            m2.boundary(Field::v1,  &chi7, BC::Neumann);

            Structure s({&m1, &m2});

            // s.inPlane11(100);
            // s.inPlane12(0);
            // s.inPlane22(100);
            s.load(1);
            s.setIterations(100);

            s.planeStrain();
            s.linear();

            s.output("plot/Data/Structure/v1", Field::v1);
            s.output("plot/Data/Structure/v2", Field::v2);
            s.output("plot/Data/Structure/z",  Field::z);

            break;
        }
        case 2: // Platte mit Loch 1
        {
            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;

            Point p0(0, 0);
            Point p1(0, 1);
            Point p2(1, 0);
            Point p3(sqrt(2)/2, sqrt(2)/2);
            Point p4(0, 2);
            Point p5(2, 0);
            Point p6(2, 2);

            Lagrange::CurveInterpolant chi1(p1, p3, p0, 1, n1);
            Lagrange::CurveInterpolant chi2(p3, p2, p0, 1, n2);
            Lagrange::CurveInterpolant chi3(p1, p4, n3);
            Lagrange::CurveInterpolant chi4(p2, p5, n3);
            Lagrange::CurveInterpolant chi5(p4, p6, n1);
            Lagrange::CurveInterpolant chi6(p5, p6, n2);
            Lagrange::CurveInterpolant chi7(p3, p6, n3);

            Membrane m1({&chi1, &chi7, &chi5, &chi3});
            Membrane m2({&chi2, &chi4, &chi6, &chi7});

            m1.boundary(Field::z,   &chi1, BC::Neumann);
            m1.boundary(Field::n12, &chi1, BC::Dirichlet);
            m1.boundary(Field::n22, &chi1, BC::Dirichlet);

            m1.boundary(Field::z,   &chi3, BC::Neumann);
            m1.boundary(Field::v2,  &chi3, BC::Neumann);
            m1.boundary(Field::v1,  &chi3, BC::Dirichlet);

            m1.boundary(Field::z,   &chi5, BC::Dirichlet);
            m1.boundary(Field::n22, &chi5, BC::Dirichlet, 1);
            m1.boundary(Field::n12, &chi5, BC::Dirichlet);

            m2.boundary(Field::z,   &chi2, BC::Neumann);
            m2.boundary(Field::n12, &chi2, BC::Dirichlet);
            m2.boundary(Field::n22, &chi2, BC::Dirichlet);
            
            m2.boundary(Field::z,   &chi4, BC::Neumann);
            m2.boundary(Field::v2,  &chi4, BC::Neumann);
            m2.boundary(Field::v1,  &chi4, BC::Dirichlet);

            m2.boundary(Field::z,   &chi6, BC::Dirichlet);
            m2.boundary(Field::n22, &chi6, BC::Dirichlet, 1);
            m2.boundary(Field::n12, &chi6, BC::Dirichlet);

            Structure s({&m1, &m2});

            s.load(1);
            // s.inPlane11(1);
            // s.inPlane12(0);
            // s.inPlane22(1);
            s.setIterations(1000);

            s.planeStrain();
            s.linear();

            s.output("plot/Data/Structure/z",  Field::z);
            s.output("plot/Data/Structure/v1", Field::v1);
            s.output("plot/Data/Structure/v2", Field::v2);
            s.output("plot/Data/Structure/n11", Field::n11);
            s.output("plot/Data/Structure/n12", Field::n12);
            s.output("plot/Data/Structure/n22", Field::n22);
            s.principalStresses("plot/Data/Structure/kartesianStresses", "plot/Data/Structure/principalStresses");
            s.principalStrains("plot/Data/Structure/kartesianDeformations", "plot/Data/Structure/kartesianStrains", "plot/Data/Structure/principalStrains");

            break;
        }
        case 3: // Platte it Loch 2
        {
            double D = 0.5; // Height of the plate in m
            double d = 0.1; // Diameter of the hole in m
            double h = 0.15; // Distance between hole and border in m (redundant)
            double t = 0.01; // Thickness of the plate in m
            double delta = 0.5*(D-d)/2; // Location of refined mesh

            double E = 2.1e11; // Youngs modulus
            double nu = 0.3; // Poissons ratio

            double P = 20000; // Load in N/m (Force divided by D)

            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;
            size_t n4 = 10;

            Point p0(0, 0);
            Point p1(0, d/2);
            Point p2(d/2, 0);
            Point p3(sqrt(2)/2*d/2, sqrt(2)/2*d/2);
            Point p4(0, D/2);
            Point p5(D/2, 0);
            Point p6(D/2, D/2);
            Point p7(d/2+delta, 0);
            Point p8(0, d/2+delta);
            Point p9((d/2+delta)*sqrt(2)/2, (d/2+delta)*sqrt(2)/2);

            Lagrange::CurveInterpolant chi1 (p1, p3, p0, d/2, n1);
            Lagrange::CurveInterpolant chi2 (p3, p2, p0, d/2, n2);
            Lagrange::CurveInterpolant chi3 (p1, p8, n3);
            Lagrange::CurveInterpolant chi4 (p2, p7, n3);
            Lagrange::CurveInterpolant chi5 (p4, p6, n1);
            Lagrange::CurveInterpolant chi6 (p5, p6, n2);
            Lagrange::CurveInterpolant chi7 (p3, p9, n3);
            Lagrange::CurveInterpolant chi8 (p9, p6, n3);
            Lagrange::CurveInterpolant chi9 (p8, p4, n4);
            Lagrange::CurveInterpolant chi10(p7, p5, n4);
            Lagrange::CurveInterpolant chi11(p8, p9, p0, d/2+delta, n4);
            Lagrange::CurveInterpolant chi12(p9, p7, p0, d/2+delta, n4);

            Membrane m1({&chi1,  &chi7,  &chi11, &chi3});
            Membrane m2({&chi2,  &chi4,  &chi12, &chi7});
            Membrane m3({&chi11, &chi8,  &chi5,  &chi9});
            Membrane m4({&chi12, &chi10, &chi6,  &chi8});

            m1.boundary(Field::n12, &chi1,  BC::Dirichlet);
            m1.boundary(Field::n22, &chi1,  BC::Dirichlet);
            m1.boundary(Field::v2,  &chi3,  BC::Neumann);
            m1.boundary(Field::v1,  &chi3,  BC::Dirichlet);

            m2.boundary(Field::n12, &chi2,  BC::Dirichlet);
            m2.boundary(Field::n22, &chi2,  BC::Dirichlet);
            m2.boundary(Field::v2,  &chi4,  BC::Neumann);
            m2.boundary(Field::v1,  &chi4,  BC::Dirichlet);

            m3.boundary(Field::v2,  &chi9,  BC::Neumann);
            m3.boundary(Field::v1,  &chi9,  BC::Dirichlet);
            m3.boundary(Field::n12, &chi5,  BC::Dirichlet);
            m3.boundary(Field::n22, &chi5,  BC::Dirichlet);

            m4.boundary(Field::v2,  &chi10, BC::Neumann);
            m4.boundary(Field::v1,  &chi10, BC::Dirichlet);
            m4.boundary(Field::n12, &chi6,  BC::Dirichlet, P);
            m4.boundary(Field::n22, &chi6,  BC::Dirichlet);

            Structure s({&m1, &m2, &m3, &m4});

            s.youngsModulus(E*t);
            s.poissonsRatio(nu);
            s.setIterations(500);
            s.planeStrain();

            s.output("plot/Data/Structure/v1",  Field::v1);
            s.output("plot/Data/Structure/v2",  Field::v2);
            s.output("plot/Data/Structure/n11", Field::n11);
            s.output("plot/Data/Structure/n12", Field::n12);
            s.output("plot/Data/Structure/n22", Field::n22);
            s.principalStresses("plot/Data/Structure/kartesianStresses", "plot/Data/Structure/principalStresses");
            s.principalStrains("plot/Data/Structure/kartesianDeformations", "plot/Data/Structure/kartesianStrains", "plot/Data/Structure/principalStrains");
            break;
        }
        case 4: // Quadrat
        {
            size_t n1 = 10;
            size_t n2 = 11;
            size_t n3 = 12;
            size_t n4 = 13;

            Point p1(0, 0);
            Point p2(0, 2);
            Point p3(0, 4);
            Point p4(2, 0);
            Point p5(2, 2);
            Point p6(2, 4);
            Point p7(4, 0);
            Point p8(4, 2);
            Point p9(4, 4);

            Lagrange::CurveInterpolant chi1 (p1, p4, n1);
            Lagrange::CurveInterpolant chi2 (p4, p5, n2);
            Lagrange::CurveInterpolant chi3 (p2, p5, n1);
            Lagrange::CurveInterpolant chi4 (p1, p2, n2);
            Lagrange::CurveInterpolant chi5 (p5, p6, n3);
            Lagrange::CurveInterpolant chi6 (p3, p6, n1);
            Lagrange::CurveInterpolant chi7 (p2, p3, n3);
            Lagrange::CurveInterpolant chi8 (p4, p7, n4);
            Lagrange::CurveInterpolant chi9 (p7, p8, n2);
            Lagrange::CurveInterpolant chi10(p5, p8, n4);
            Lagrange::CurveInterpolant chi11(p8, p9, n3);
            Lagrange::CurveInterpolant chi12(p6, p9, n4);

            Membrane m1({&chi1,  &chi2,  &chi3,  &chi4});
            Membrane m2({&chi3,  &chi5,  &chi6,  &chi7});
            Membrane m3({&chi8,  &chi9,  &chi10, &chi2});
            Membrane m4({&chi10, &chi11, &chi12, &chi5});

            m1.boundary(Field::v2,  &chi1,  BC::Dirichlet);
            m1.boundary(Field::v1,  &chi1,  BC::Neumann);
            m1.boundary(Field::z,   &chi1,  BC::Neumann);
            m1.boundary(Field::v2,  &chi4,  BC::Neumann);
            m1.boundary(Field::v1,  &chi4,  BC::Dirichlet);
            m1.boundary(Field::z,   &chi4,  BC::Neumann);

            m2.boundary(Field::v1,  &chi7,  BC::Dirichlet);
            m2.boundary(Field::v2,  &chi7,  BC::Neumann);
            m2.boundary(Field::z,   &chi7,  BC::Neumann);
            m2.boundary(Field::n12, &chi6,  BC::Dirichlet);
            m2.boundary(Field::n22, &chi6,  BC::Dirichlet, 1);
            m2.boundary(Field::z,   &chi6,  BC::Dirichlet);

            m3.boundary(Field::n11, &chi9,  BC::Dirichlet, 1);
            m3.boundary(Field::n12, &chi9,  BC::Dirichlet);
            m3.boundary(Field::z,   &chi9,  BC::Dirichlet);
            m3.boundary(Field::v1,  &chi8,  BC::Neumann);
            m3.boundary(Field::v2,  &chi8,  BC::Dirichlet);
            m3.boundary(Field::z,   &chi8,  BC::Neumann);

            m4.boundary(Field::n12, &chi11, BC::Dirichlet);
            m4.boundary(Field::n11, &chi11, BC::Dirichlet, 1);
            m4.boundary(Field::z,   &chi11, BC::Dirichlet);
            m4.boundary(Field::n22, &chi12, BC::Dirichlet, 1);
            m4.boundary(Field::n12, &chi12, BC::Dirichlet);
            m4.boundary(Field::z,   &chi12, BC::Dirichlet);

            Structure s({&m4, &m3, &m2, &m1});

            s.load(1);

            // s.inPlane11(1);
            // s.inPlane12(0);
            // s.inPlane22(1);

            s.setIterations(1000);
            s.planeStrain();
            s.linear();

            s.output("plot/Data/Structure/v1", Field::v1);
            s.output("plot/Data/Structure/v2", Field::v2);
            s.output("plot/Data/Structure/z",  Field::z);
            s.principalStresses("plot/Data/Structure/kartesianStresses", "plot/Data/Structure/principalStresses");
            s.principalStrains("plot/Data/Structure/kartesianDeformations", "plot/Data/Strucutre/kartesianStrains", "plot/Data/Structure/principalStrains");

            break;
        }
        case 5: // L
        {
            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;
            size_t n4 = 10;

            Point p1(0, 0);
            Point p2(0, 2);
            Point p3(0, 4);
            Point p4(2, 0);
            Point p5(2, 2);
            Point p6(2, 4);
            Point p7(4, 0);
            Point p8(4, 2);

            Lagrange::CurveInterpolant chi1 (p1, p4, n1);
            Lagrange::CurveInterpolant chi2 (p4, p5, n2);
            Lagrange::CurveInterpolant chi3 (p2, p5, n1);
            Lagrange::CurveInterpolant chi4 (p1, p2, n2);
            Lagrange::CurveInterpolant chi5 (p5, p6, n3);
            Lagrange::CurveInterpolant chi6 (p3, p6, n1);
            Lagrange::CurveInterpolant chi7 (p2, p3, n3);
            Lagrange::CurveInterpolant chi8 (p4, p7, n4);
            Lagrange::CurveInterpolant chi9 (p7, p8, n2);
            Lagrange::CurveInterpolant chi10(p5, p8, n4);

            Membrane m1({&chi1, &chi2, &chi3,  &chi4});
            Membrane m2({&chi3, &chi5, &chi6,  &chi7});
            Membrane m3({&chi8, &chi9, &chi10, &chi2});

            m1.boundary(Field::v2,  &chi1,  BC::Dirichlet);
            m1.boundary(Field::v1,  &chi1,  BC::Neumann);
            m1.boundary(Field::z,   &chi1,  BC::Neumann);
            m1.boundary(Field::v2,  &chi4,  BC::Neumann);
            m1.boundary(Field::v1,  &chi4,  BC::Dirichlet);
            m1.boundary(Field::z,   &chi4,  BC::Neumann);

            m2.boundary(Field::v1,  &chi7,  BC::Dirichlet);
            m2.boundary(Field::v2,  &chi7,  BC::Neumann);
            m2.boundary(Field::z,   &chi7,  BC::Neumann);
            m2.boundary(Field::n12, &chi6,  BC::Dirichlet);
            m2.boundary(Field::n22, &chi6,  BC::Dirichlet, 1);
            m2.boundary(Field::z,   &chi6,  BC::Dirichlet);
            m2.boundary(Field::n11, &chi5,  BC::Dirichlet, 1);
            m2.boundary(Field::n12, &chi5,  BC::Dirichlet);
            m2.boundary(Field::z,   &chi5,  BC::Dirichlet);

            m3.boundary(Field::n12, &chi10, BC::Dirichlet);
            m3.boundary(Field::n22, &chi10, BC::Dirichlet, 1);
            m3.boundary(Field::z,   &chi10, BC::Dirichlet);
            m3.boundary(Field::n11, &chi9,  BC::Dirichlet, 1);
            m3.boundary(Field::n12, &chi9,  BC::Dirichlet);
            m3.boundary(Field::z,   &chi9,  BC::Dirichlet);
            m3.boundary(Field::v1,  &chi8,  BC::Neumann);
            m3.boundary(Field::v2,  &chi8,  BC::Dirichlet);
            m3.boundary(Field::z,   &chi8,  BC::Neumann);

            Structure s({&m1, &m2, &m3});

            // s.inPlane11(100);
            // s.inPlane12(0);
            // s.inPlane22(100);
            s.load(1);
            s.setIterations(1000);
            s.planeStrain();
            s.linear();

            s.output("plot/Data/Structure/v1", Field::v1);
            s.output("plot/Data/Structure/v2", Field::v2);
            s.output("plot/Data/Structure/z",  Field::z);

            break;
        }
        case 6: // Kreis
        {
            size_t n1 = 10;
            size_t n2 = 10;
            size_t n3 = 10;

            double r = 1;

            Point p1(0, 0);
            Point p2(r/2, 0);
            Point p3(0.45*r, 0.45*r);
            Point p4(0, r/2);
            Point p5(r, 0);
            Point p6(0, r);
            Point p7(r/sqrt(2), r/sqrt(2));

            Lagrange::CurveInterpolant chi1(p1, p2, n1);
            Lagrange::CurveInterpolant chi2(p2, p3, n2);
            Lagrange::CurveInterpolant chi3(p4, p3, n1);
            Lagrange::CurveInterpolant chi4(p1, p4, n2);
            Lagrange::CurveInterpolant chi5(p2, p5, n3);
            Lagrange::CurveInterpolant chi6(p4, p6, n3);
            Lagrange::CurveInterpolant chi7(p3, p7, n3);
            Lagrange::CurveInterpolant chi8(p5, p7, p1, r, n2);
            Lagrange::CurveInterpolant chi9(p6, p7, p1, r, n1);

            Membrane m1({&chi1, &chi2, &chi3, &chi4});
            Membrane m2({&chi5, &chi8, &chi7, &chi2});
            Membrane m3({&chi3, &chi7, &chi9, &chi6});

            m1.boundary(Field::z,   &chi1, BC::Neumann);
            m1.boundary(Field::z,   &chi4, BC::Neumann);
            m2.boundary(Field::z,   &chi5, BC::Neumann);
            m2.boundary(Field::z,   &chi8, BC::Dirichlet);
            m3.boundary(Field::z,   &chi6, BC::Neumann);
            m3.boundary(Field::z,   &chi9, BC::Dirichlet);
 
            m1.boundary(Field::v1,  &chi1, BC::Neumann);
            m1.boundary(Field::v1,  &chi4, BC::Dirichlet);
            m2.boundary(Field::v1,  &chi5, BC::Neumann);
            m3.boundary(Field::v1,  &chi6, BC::Dirichlet);
 
            m1.boundary(Field::v2,  &chi1, BC::Dirichlet);
            m1.boundary(Field::v2,  &chi4, BC::Neumann);
            m2.boundary(Field::v2,  &chi5, BC::Dirichlet);
            m3.boundary(Field::v2,  &chi6, BC::Neumann);

            m2.boundary(Field::n11, &chi8, BC::Dirichlet, 10);
            m2.boundary(Field::n12, &chi8, BC::Dirichlet);
            m3.boundary(Field::n12, &chi9, BC::Dirichlet);
            m3.boundary(Field::n22, &chi9, BC::Dirichlet, 10);

            Structure s({&m1, &m2, &m3});

            s.load(1);
            s.setIterations(200);
            s.planeStrain();
            s.linear();
            s.semilinear();
            s.nonlinear();

            s.output("plot/Data/Structure/v1",  Field::v1);
            s.output("plot/Data/Structure/v2",  Field::v2);
            s.output("plot/Data/Structure/z",   Field::z);
            s.output("plot/Data/Structure/n11", Field::n11);
            s.output("plot/Data/Structure/n12", Field::n12);
            s.output("plot/Data/Structure/n22", Field::n22);
            s.principalStresses("plot/Data/Structure/kartesianStresses", "plot/Data/Structure/principalStresses");
            s.principalStrains("plot/Data/Structure/kartesianDeformations", "plot/Data/Structure/kartesianStrains", "plot/Data/Structure/principalStrains");

            break;
        }
    }
}