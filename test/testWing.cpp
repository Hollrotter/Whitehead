#include "Wing.hpp"

int main()
{
    switch (2)
    {
        case 0: // Rectangle
        {
            size_t nx = 20;
            size_t ny = 50;

            double l = 2;
            double b = 10;

            Point p1(0,-b/2);
            Point p2(l,-b/2);
            Point p3(l, b/2);
            Point p4(0, b/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p3, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, nx);
            Lagrange::CurveInterpolant chi4(p4, p1, ny);

            Wing w({&chi1, &chi2, &chi3, &chi4});

            w.pitch(2);

            w.boundary(&chi1, BC::Dirichlet);
            w.boundary(&chi2, BC::Neumann);
            w.boundary(&chi3, BC::Dirichlet);
            w.boundary(&chi4, BC::Dirichlet);

            w.linear();

            double area = w.get_area();
            arma::vec cL = w.get_lift()   / area;
            arma::vec cM = w.get_moment() / area/l;

            std::cout << "A  = " << area   << '\n';
            std::cout << "cL = " << cL.t() << '\n';
            std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/flat");
            break;
        }
        case 1: // Symmetry
        {
            /**
             * Basic test for the symmetry in the y-direction.
             * The result must be the same like in the
             * rectangle test case.
             * The results do not fit exactly, because the
             * two cases have different grids. For a high
             * number of nodes, we get exact agreement.
             * For nx = 20 and ny = 25, we got
             * cL = 0.1677 and cM =-0.0424
             * and the results do agree.
            */

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
            w(Symmetry::y);

            w.boundary(&chi1, BC::Neumann);
            w.boundary(&chi2, BC::Neumann);
            w.boundary(&chi3, BC::Dirichlet);
            w.boundary(&chi4, BC::Dirichlet);

            w.linear();

            double area = w.get_area();
            arma::vec cL = w.get_lift()   / area;
            arma::vec cM = w.get_moment() / area/l;

            std::cout << "A  = " << area   << '\n';
            std::cout << "cL = " << cL.t() << '\n';
            std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/symmetry");
            break;
        }
        case 2: // Rotated square
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

            size_t n = 20;
            size_t m = 10;

            Point p1(-1, 0);
            Point p2( 1, 0);
            Point p3( 1, 1);
            Point p4(-1, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, m);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, m);

            arma::mat z1 = 0.1*cos(arma::datum::pi/4*(1-Chebyshev::gaussLobatto(n))) * cos(arma::datum::pi/4*(1+Chebyshev::gaussLobatto(m))).t();
            arma::mat z2 = 0.1*cos(arma::datum::pi/4*(1+Chebyshev::gaussLobatto(m))) * cos(arma::datum::pi/4*(1+Chebyshev::gaussLobatto(n))).t();
            arma::mat z3 = arma::flipud(arma::fliplr(z1));
            arma::mat z4 = arma::flipud(arma::fliplr(z2));

            std::array<Wing, 4> w({Wing(z1, {&chi1, &chi2, &chi3, &chi4}),
                                   Wing(z2, {&chi2, &chi3, &chi4, &chi1}),
                                   Wing(z3, {&chi3, &chi4, &chi1, &chi2}),
                                   Wing(z4, {&chi4, &chi1, &chi2, &chi3})});
            Wake wk(&chi2);
            for (size_t i = 0; i < 4; i++)
            {
                w[i].pitch(3);
                w[i].boundary(&chi1, BC::Neumann);
                w[i].boundary(&chi2, BC::Neumann);
                w[i].boundary(&chi3, BC::Dirichlet);
                w[i].boundary(&chi4, BC::Dirichlet);
                w[i].wake(&wk);
                w[i](Symmetry::y);
                w[i].nonlinear();
                w[i].output("plot/Data/Wing/square"+std::to_string(i));
            }

            break;
        }
        case 3: // Elliptical Wing
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

            w.boundary(&chi1, BC::Neumann);
            w.boundary(&chi2, BC::Derivative_x);
            w.boundary(&chi3, BC::Dirichlet);
            w.boundary(&chi4, BC::Dirichlet);

            w.linear();

            double area = w.get_area();
            arma::vec cL = w.get_lift()   / area;
            arma::vec cM = w.get_moment() / area/c;

            std::cout << "A  = " << area   << '\n';
            std::cout << "cL = " << cL.t() << '\n';
            std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/EllipticalWing");
            break;
        }
        case 4: // Nonlinear
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

            double area = 2*w.get_area();
            arma::vec cL = w.get_lift()   / area;
            arma::vec cM = w.get_moment() / area/l;

            std::cout << "A  = " << area   << '\n';
            std::cout << "cL = " << cL.t() << '\n';
            std::cout << "cM = " << cM.t() << '\n';

            w.output("plot/Data/Wing/nonlinear");
            break;
        }
    }
}