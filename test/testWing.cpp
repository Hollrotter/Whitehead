#include "Wing.hpp"

int main()
{
    switch (0)
    {
        case 0: // Rectangle
        {
            size_t nx = 7;
            size_t ny = 7;

            double l = 2;
            double b = 10;

            Point p1(0, 0);
            Point p2(l, 0);
            Point p3(l, b/2);
            Point p4(0, b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);

            Wing w({&chi1, &chi2, &chi3, &chi4});

            w.pitch(2);

            w.boundary(&chi1, BC::Neumann);
            w.boundary(&chi2, BC::Neumann);
            w.boundary(&chi3, BC::Dirichlet);
            w.boundary(&chi4, BC::Dirichlet);

            w.linear();

            w.output("plot/Data/Wing/flat");
            break;
        }
    }
}