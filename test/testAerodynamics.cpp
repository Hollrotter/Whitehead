#include "Aerodynamics.hpp"

int main()
{
    switch (0)
    {
        case 0: // Rectangle (divided at y=0)
        {
            size_t nx = 7;
            size_t ny = 7;

            double l = 2;
            double b = 4;

            Point p1(-l/2, 0);
            Point p2( l/2, 0);
            Point p3( l/2, b/2);
            Point p4(-l/2, b/2);
            Point p5(-l/2,-b/2);
            Point p6( l/2,-b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);
            Lagrange::CurveInterpolant chi5(p1, p5, ny);
            Lagrange::CurveInterpolant chi6(p2, p6, ny);
            Lagrange::CurveInterpolant chi7(p5, p6, nx);

            Wing w1({&chi1, &chi2, &chi3, &chi4});
            Wing w2({&chi7, &chi6, &chi1, &chi5});

            Aerodynamics a({&w1, &w2});
            a.pitch(2);

            w1.boundary(&chi2, BC::Neumann);
            w1.boundary(&chi3, BC::Dirichlet);
            w1.boundary(&chi4, BC::Dirichlet);

            w2.boundary(&chi7, BC::Dirichlet);
            w2.boundary(&chi6, BC::Neumann);
            w2.boundary(&chi5, BC::Dirichlet);

            a.setIterations(100);
            a.linear();

            a.output("plot/Data/Aerodynamics/flat");
            break;
        }
    }
}