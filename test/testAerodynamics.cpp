#include "Aerodynamics.hpp"

int main()
{
    switch (3)
    {
        case 0: // Rectangle (divided at y=0)
        {
            size_t nx = 25;
            size_t ny = 20;

            double AR = 5;
            double l = 2;
            double b = l * AR;

            double alpha = 5;

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
            a.pitch(alpha);

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
            double cLalpha = cL/(alpha*arma::datum::pi/180);

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

            a.boundary(&chi2, BC::Kutta);
            a.boundary(&chi3, BC::Dirichlet);
            a.boundary(&chi4, BC::Dirichlet);
            a.boundary(&chi5, BC::Dirichlet);
            a.boundary(&chi6, BC::Kutta);
            a.boundary(&chi7, BC::Dirichlet);

            a.setIterations(100);
            a.nonlinear();

            a.output("plot/Data/Aerodynamics/nonsmooth");
            break;
        }
        case 3: // Panel
        {
            /**
             * Rectangular Wing split up into an arbitrary number of Panels.
             * This will in practice never be done. The purpose of this test
             * is to see the results with multiple surfaces and interfaces.
             */
            double x_min =-1;
            double x_max = 1;
            double y_min =-1;
            double y_max = 1;

            const size_t nx = 6; // Number of Panels in x-direction
            const size_t ny = 6; // Number of Panels in y-direction

            double dx = (x_max-x_min)/nx;
            double dy = (y_max-y_min)/ny;

            size_t n1 = 10; // Number of nodes per panel in x-direction
            size_t n2 = 10; // Number of nodes per panel in y-direction

            std::array<std::array<Point, ny+1>, nx+1> points;
            std::array<std::array<Lagrange::CurveInterpolant, ny+1>, nx> chiH;
            std::array<std::array<Lagrange::CurveInterpolant, ny>, nx+1> chiV;
            std::vector<Wing> wings;
            wings.reserve(nx*ny);

            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    points[i][j]     = Point(x_min+i*dx, y_min+j*dy);
                    points[i+1][j]   = Point(x_min+(i+1)*dx, y_min+j*dy);
                    points[i][j+1]   = Point(x_min+i*dx, y_min+(j+1)*dy);
                    points[i+1][j+1] = Point(x_min+(i+1)*dx, y_min+(j+1)*dy);
                }
            
            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    chiH[i][j]   = Lagrange::CurveInterpolant(points[i][j],   points[i+1][j],   n1);
                    chiH[i][j+1] = Lagrange::CurveInterpolant(points[i][j+1], points[i+1][j+1], n1);
                    chiV[i][j]   = Lagrange::CurveInterpolant(points[i][j],   points[i][j+1],   n2);
                    chiV[i+1][j] = Lagrange::CurveInterpolant(points[i+1][j], points[i+1][j+1], n2);
                    Wing w({&chiH[i][j], &chiV[i+1][j], &chiH[i][j+1], &chiV[i][j]});
                    wings.emplace_back(w);
                }

            std::vector<Wing*> wings_ptr;
            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                    wings_ptr.push_back(&wings[i+j*nx]);

            Aerodynamics panel(wings_ptr);
            panel.setIterations(500);
            panel.setlambda(10);
            panel.pitch(5);

            for (size_t j = 0; j < ny; j++)
            {
                panel.boundary(&chiV[0][j],  BC::Dirichlet);
                panel.boundary(&chiV[nx][j], BC::Dirichlet);
            }

            for (size_t i = 0; i < nx; i++)
            {
                panel.boundary(&chiH[i][0],  BC::Dirichlet);
                panel.boundary(&chiH[i][ny], BC::Dirichlet);
            }

            panel.linear();

            panel.output("plot/Data/Aerodynamics/panel");
            break;
        }
    }
}