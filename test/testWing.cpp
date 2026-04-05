#include "Wing.hpp"

int main()
{
    switch (3)
    {
        case 0: // Rectangle
        {
            size_t nx = 15;
            size_t ny = 25;

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

            double area = 2*w.get_area();
            arma::vec cL = w.get_lift()   / area;
            arma::vec cM = w.get_moment() / area/l;

            std::cout << "A  = " << area   << '\n';
            std::cout << "cL = " << cL.t() << '\n';
            std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/flat");
            break;
        }
        case 1: // Rotated square
        {
            /**
             * In this test case, we solve the same proble of a
             * square wing four times, but rotated by a 90 degree
             * angle. In this basic test, the computational grid is
             * identical with the phyisical grid, except the rotation.
             * The result must be the same for all four cases. The main
             * takeaway is, that the transformation is in principle doing
             * what it is supposed to do.
             */

            size_t n = 15;

            Point p1(-1,-1);
            Point p2( 1,-1);
            Point p3( 1, 1);
            Point p4(-1, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            std::array<Wing, 4> w({Wing({&chi1, &chi2, &chi3, &chi4}),
                                   Wing({&chi2, &chi3, &chi4, &chi1}),
                                   Wing({&chi3, &chi4, &chi1, &chi2}),
                                   Wing({&chi4, &chi1, &chi2, &chi3})});

            for (size_t i = 0; i < 4; i++)
            {
                w[i].pitch(2);
                w[i].boundary(&chi1, BC::Dirichlet);
                w[i].boundary(&chi2, BC::Neumann);
                w[i].boundary(&chi3, BC::Dirichlet);
                w[i].boundary(&chi4, BC::Dirichlet);
                w[i].linear();
                w[i].output("plot/Data/Wing/square"+std::to_string(i));
            }

            break;
        }
        case 2: // Elliptical Wing
        {
            double b = 10;
            double c = 1.39;

            size_t nx = 10;
            size_t ny = 15;

            arma::vec x1 = c/2*(1+Chebyshev::gaussLobatto(nx));
            arma::vec x2 = 0.9*b/4*(1 + Chebyshev::gaussLobatto(ny));

            arma::mat y = arma::ones(x1.size(), 1)*x2.t();
            arma::mat x = (x1-c/4)*sqrt(1-pow(x2/(b/2), 2)).t();

            Lagrange::CurveInterpolant chi1(x.col(0),        y.col(0));
            Lagrange::CurveInterpolant chi2(x.row(nx-1).t(), y.row(nx-1).t());
            Lagrange::CurveInterpolant chi3(x.col(ny-1),     y.col(ny-1));
            Lagrange::CurveInterpolant chi4(x.row(0).t(),    y.row(0).t());

            Wing w({&chi1, &chi2, &chi3, &chi4});

            w.pitch(5);

            w.boundary(Direction::S, BC::Neumann);
            w.boundary(Direction::N, BC::Dirichlet);
            w.boundary(Direction::W, BC::Dirichlet);
            w.boundary(Direction::E, BC::Neumann);

            w.linear();

            w.output("plot/Data/Wing/EllipticalWing");
            break;
        }
        case 3: // Nonlinear
        {
            size_t nx = 10;
            size_t ny = 10;

            double l = 1;
            double b = 10;

            Point p1(-l/4, 0);
            Point p2(3*l/4, 0);
            Point p3(3*l/4, b/2);
            Point p4(-l/4, b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);

            arma::mat z = 0.2 * cos(arma::datum::pi/2* Chebyshev::gaussLobatto(nx))
                               * cos(arma::datum::pi/4*(Chebyshev::gaussLobatto(ny)+1)).t();
            Wing w(z, {&chi1, &chi2, &chi3, &chi4});

            w.pitch(5);

            w.boundary(&chi1, BC::Neumann);
            w.boundary(&chi2, BC::Neumann);
            w.boundary(&chi3, BC::Dirichlet);
            w.boundary(&chi4, BC::Dirichlet);

            w.nonlinear();

            // double area = 2*w.get_area();
            // arma::vec cL = w.get_lift()   / area;
            // arma::vec cM = w.get_moment() / area/l;

            // std::cout << "A  = " << area   << '\n';
            // std::cout << "cL = " << cL.t() << '\n';
            // std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/nonlinear");
            break;
        }
    }
}