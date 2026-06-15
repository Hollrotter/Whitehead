#include "Aerodynamics.hpp"

int main()
{
    switch (2)
    {
        case 0: // Rectangle (divided at y=0)
        {
            size_t nx = 10;
            size_t ny = 10;

            double l = 2;
            double AR = 5;
            double b = l * AR;

            Point p1(0, 0);
            Point p2(l, 0);
            Point p3(l, b/2);
            Point p4(0, b/2);
            Point p5(0,-b/2);
            Point p6(l,-b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);
            Lagrange::CurveInterpolant chi5(p1, p5, ny);
            Lagrange::CurveInterpolant chi6(p2, p6, ny);
            Lagrange::CurveInterpolant chi7(p5, p6, nx);

            Wing w1({&chi1, &chi2, &chi3, &chi4});
            Wing w2({&chi7, &chi6, &chi1, &chi5});

            Wake wake1(&chi2);
            Wake wake2(&chi6);

            Aerodynamics a({&w1, &w2});
            a.pitch(5);

            a.boundary(&chi2, BC::Neumann);
            a.boundary(&chi3, BC::Dirichlet);
            a.boundary(&chi4, BC::Dirichlet);
            a.boundary(&chi5, BC::Dirichlet);
            a.boundary(&chi6, BC::Neumann);
            a.boundary(&chi7, BC::Dirichlet);

            a.wake(&wake1);
            a.wake(&wake2);

            a.setIterations(100);
            a.linear();

            double area = a.get_area();
            double cL   = a.get_lift()   / area;
            double cM   = a.get_moment() / area/l;

            std::cout << "A  = " << area << '\n';
            std::cout << "cL = " << cL   << '\n';
            std::cout << "cM = " << cM   << '\n';

            double a0 = arma::datum::tau;
            double cLalpha_A1 = a0/(1 + a0/arma::datum::pi/AR);
            double cLalpha_A2 = a0/(1 + 1.024*a0/arma::datum::pi/AR);
            double cLalpha_A3 = a0/(sqrt(1 + pow(a0/arma::datum::pi/AR, 2)) + a0/arma::datum::pi/AR);
            double cLalpha = cL/(5*arma::datum::pi/180);

            std::cout << "cLalpha    = " << cLalpha   << '\n';
            std::cout << "cLalpha_A1 = " << cLalpha_A1 << '\n';
            std::cout << "cLalpha_A2 = " << cLalpha_A2 << '\n';
            std::cout << "cLalpha_A3 = " << cLalpha_A3 << '\n';

            a.output("plot/Data/Aerodynamics/flat");
            break;
        }
        case 1: // More complex shape
        {
            Point p1(-1, 0);
            Point p2( 1, 0);
            Point p3( 0, 1);
            Point p4( 0,-1);
            Point p5(0.5, 4);
            Point p6(1, 4);
            Point p7(0.5,-4);
            Point p8(1,-4);

            size_t n = 15;

            Lagrange::CurveInterpolant chi1 (p1, p3, n);
            Lagrange::CurveInterpolant chi2 (p2, p3, n);
            Lagrange::CurveInterpolant chi3 (p2, p4, n);
            Lagrange::CurveInterpolant chi4 (p1, p4, n);
            Lagrange::CurveInterpolant chi5 (p3, p5, n);
            Lagrange::CurveInterpolant chi6 (p5, p6, n);
            Lagrange::CurveInterpolant chi7 (p2, p6, n);
            Lagrange::CurveInterpolant chi8 (p4, p7, n);
            Lagrange::CurveInterpolant chi9 (p2, p8, n);
            Lagrange::CurveInterpolant chi10(p7, p8, n);

            Wing w1({&chi1, &chi2, &chi3,  &chi4});
            Wing w2({&chi5, &chi6, &chi7,  &chi2});
            Wing w3({&chi3, &chi9, &chi10, &chi8});

            Aerodynamics a({&w1, &w2, &w3});

            a.boundary(&chi1,  BC::Dirichlet);
            a.boundary(&chi4,  BC::Dirichlet);
            a.boundary(&chi5,  BC::Dirichlet);
            a.boundary(&chi8,  BC::Dirichlet);
            a.boundary(&chi6,  BC::Dirichlet);
            a.boundary(&chi7,  BC::Neumann);
            a.boundary(&chi9,  BC::Neumann);
            a.boundary(&chi10, BC::Dirichlet);

            a.pitch(1);
            a.setIterations(75);
            a.linear();

            a.output("plot/Data/Aerodynamics/v");
            break;
        }
        case 2: // Nonsmooth Surface
        {
            size_t nx = 10;
            size_t ny = 15;

            double l = 2;
            double b = 4;

            Point p1(0, 0);
            Point p2(l, 0);
            Point p3(l, b/2);
            Point p4(0, b/2);
            Point p5(0,-b/2);
            Point p6(l,-b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);
            Lagrange::CurveInterpolant chi5(p1, p5, ny);
            Lagrange::CurveInterpolant chi6(p2, p6, ny);
            Lagrange::CurveInterpolant chi7(p5, p6, nx);

            arma::mat z1 = cos(arma::datum::pi/2*Chebyshev::gaussLobatto(nx))*cos(arma::datum::pi/2*Chebyshev::gaussLobatto(ny)).t();
            arma::mat z2 = cos(arma::datum::pi/2*Chebyshev::gaussLobatto(nx))*cos(arma::datum::pi/2*Chebyshev::gaussLobatto(ny)).t();

            Wing w1(z1, {&chi1, &chi2, &chi3, &chi4});
            Wing w2(z2, {&chi7, &chi6, &chi1, &chi5});

            Aerodynamics a({&w1, &w2});
            a.pitch(2);

            a.boundary(&chi2, BC::Derivative_x);
            a.boundary(&chi3, BC::Dirichlet);
            a.boundary(&chi4, BC::Dirichlet);
            a.boundary(&chi5, BC::Dirichlet);
            a.boundary(&chi6, BC::Derivative_x);
            a.boundary(&chi7, BC::Dirichlet);

            a.setIterations(100);
            a.nonlinear();

            a.output("plot/Data/Aerodynamics/nonsmooth");
            break;
        }
    }
}