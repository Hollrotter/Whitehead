#include "Membrane.hpp"

int main()
{
    switch(3)
    {
        case 0: // Parallelogramm
        {
            double alpha = 0 * arma::datum::pi/180;
            size_t nx = 20;
            size_t ny = 20;

            Point p1(0, 0);
            Point p2(2, 0);
            Point p3(2+2*tan(alpha), 2);
            Point p4(2*tan(alpha), 2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx); // Unten, South
            Lagrange::CurveInterpolant chi2(p2, p3, ny); // Rechts, East
            Lagrange::CurveInterpolant chi3(p4, p3, nx); // Oben, North
            Lagrange::CurveInterpolant chi4(p1, p4, ny); // Links, West

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            // m.inPlane11(1);
            // m.inPlane12(0);
            // m.inPlane22(1);

            m.boundary(Field::z,   &chi1, BC::Neumann);
            m.boundary(Field::z,   &chi2, BC::Dirichlet);
            m.boundary(Field::z,   &chi3, BC::Dirichlet);
            m.boundary(Field::z,   &chi4, BC::Neumann);

            m.boundary(Field::v2,  &chi1, BC::Dirichlet);
            m.boundary(Field::n12, &chi2, BC::Dirichlet);
            m.boundary(Field::n22, &chi3, BC::Dirichlet, 10);
            m.boundary(Field::v2,  &chi4, BC::Neumann);

            m.boundary(Field::v1,  &chi1, BC::Neumann);
            m.boundary(Field::n11, &chi2, BC::Dirichlet, 10);
            m.boundary(Field::n12, &chi3, BC::Dirichlet);
            m.boundary(Field::v1,  &chi4, BC::Dirichlet);

            m.planeStrain();
            m.load(1);

            m.iterations(100);
            // m.substepControl(10);
            m.nonlinear();

            m.output("plot/Data/Membrane/z",   Field::z);
            m.output("plot/Data/Membrane/v1",  Field::v1);
            m.output("plot/Data/Membrane/v2",  Field::v2);
            m.output("plot/Data/Membrane/n11", Field::n11);
            m.output("plot/Data/Membrane/n22", Field::n22);
            m.output("plot/Data/Membrane/n12", Field::n12);

            break;
        }
        case 1: // Ellipsenflügel
        {
            double b = 10;
            double c = 1.39;

            size_t nx = 20;
            size_t ny = 30;

            arma::vec x1 = c/2*(1+Chebyshev::gaussLobatto(nx));
            arma::vec x2 = 0.9*b/4*(1+Chebyshev::gaussLobatto(ny));

            arma::mat y = arma::ones(x1.size(), 1)*x2.t();
            arma::mat x = (x1-c/4)*sqrt(1-pow(x2/(b/2), 2)).t();

            // Lagrange::CurveInterpolant chi1(x.row(0).t(),    y.row(0).t());
            // Lagrange::CurveInterpolant chi2(x.col(ny-1),     y.col(ny-1));
            // Lagrange::CurveInterpolant chi3(x.row(nx-1).t(), y.row(nx-1).t());
            // Lagrange::CurveInterpolant chi4(x.col(0),        y.col(0));

            Lagrange::CurveInterpolant chi1(x.col(0),        y.col(0));
            Lagrange::CurveInterpolant chi2(x.row(nx-1).t(), y.row(nx-1).t());
            Lagrange::CurveInterpolant chi3(x.col(ny-1),     y.col(ny-1));
            Lagrange::CurveInterpolant chi4(x.row(0).t(),    y.row(0).t());

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::z,   &chi2, BC::Dirichlet);
            m.boundary(Field::z,   &chi4, BC::Dirichlet);
            m.boundary(Field::z,   &chi1, BC::Neumann);
            m.boundary(Field::z,   &chi3, BC::Dirichlet);

            m.boundary(Field::n11, &chi2, BC::Dirichlet, 100);
            m.boundary(Field::v1,  &chi4, BC::Dirichlet);
            m.boundary(Field::n12, &chi1, BC::Dirichlet);
            m.boundary(Field::n12, &chi3, BC::Dirichlet);

            m.boundary(Field::n12, &chi2, BC::Dirichlet);
            m.boundary(Field::n12, &chi4, BC::Dirichlet);
            m.boundary(Field::v2,  &chi1, BC::Dirichlet);
            m.boundary(Field::n22, &chi3, BC::Dirichlet, 100);

            //m.inPlane1(0.01*arma::ones(15, 25));
            //m.inPlane2(0.01*arma::ones(15, 25));

            m.planeStrain();
            m.load(100*cos(arma::datum::pi/2*Chebyshev::gaussLobatto(x1.size()))*cos(arma::datum::pi/4*(1+Chebyshev::gaussLobatto(x2.size()))).t());

            // m.linear();
            // m.output("plot/Data/z_linear", Field::z);

            // m.substepControl(10);
            // m.iterations(100);

            m.nonlinear();
            m.output("plot/Data/Membrane/z",   Field::z);
            m.output("plot/Data/Membrane/v1",  Field::v1);
            m.output("plot/Data/Membrane/v2",  Field::v2);
            m.output("plot/Data/Membrane/n11", Field::n11);
            m.output("plot/Data/Membrane/n22", Field::n22);
            m.output("plot/Data/Membrane/n12", Field::n12);

            break;
        }
        case 2: // Ring
        {
            size_t n1 = 30;
            size_t n2 = 30;

            Point p0(0, 0);
            Point p1(0, 1);
            Point p2(0, 2);
            Point p3(1, 0);
            Point p4(2, 0);

            Lagrange::CurveInterpolant chi1(p1, p3, p0, 1, n1);
            Lagrange::CurveInterpolant chi2(p3, p4, n2);
            Lagrange::CurveInterpolant chi3(p2, p4, p0, 2, n1);
            Lagrange::CurveInterpolant chi4(p1, p2, n2);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            // m.boundary(Field::z,   Direction::N, BC::Dirichlet);
            // m.boundary(Field::z,   Direction::S, BC::Dirichlet);
            // m.boundary(Field::z,   Direction::W, BC::Neumann);
            // m.boundary(Field::z,   Direction::E, BC::Neumann);

            m.boundary(Field::n22, &chi1, BC::Dirichlet);
            m.boundary(Field::n22, &chi3, BC::Dirichlet, 100);
            m.boundary(Field::v2,  &chi2, BC::Neumann);
            m.boundary(Field::v2,  &chi4, BC::Neumann);

            m.boundary(Field::n12, &chi1, BC::Dirichlet);
            m.boundary(Field::n12, &chi3, BC::Dirichlet);
            m.boundary(Field::v1,  &chi2, BC::Dirichlet);
            m.boundary(Field::v1,  &chi4, BC::Dirichlet);

            m.planeStrain();
            // m.load(20*cos(Chebyshev::gaussLobatto(x1.size()))*cos(Chebyshev::gaussLobatto(x2.size()).t()));
            // m.substepControl(1);
            // m.linear();

            // m.iterations(20);
            // m.semilinear();

            // m.relaxationFactor(0.005);
            // m.nonlinear();

            // m.output("plot/Data/z",   Field::z);
            m.output("plot/Data/Membrane/v1",  Field::v1);
            m.output("plot/Data/Membrane/v2",  Field::v2);
            m.output("plot/Data/Membrane/n11", Field::n11);
            m.output("plot/Data/Membrane/n22", Field::n22);
            m.output("plot/Data/Membrane/n12", Field::n12);
            m.principalStresses("plot/Data/Membrane/kartesianStresses", "plot/Data/Membrane/principalStresses");
            m.principalStrains("plot/Data/Membrane/kartesianDeformations", "plot/Data/Membrane/kartesianStrains", "plot/Data/Membrane/principalStrains");
            break;
        }
        case 3: // Kreis
        {
            size_t nx = 30;
            size_t ny = 30;

            double r = 1;

            Point p1(0, 0);
            Point p2(r, 0);
            Point p3(0, r);
            Point p4(1.05*r*sqrt(2)/2, 1.05*r*sqrt(2)/2);

            Lagrange::CurveInterpolant chi1(p1, p2, nx);
            Lagrange::CurveInterpolant chi2(p2, p4, p1, 1.05*r, ny);
            Lagrange::CurveInterpolant chi3(p3, p4, p1, 1.05*r, nx);
            Lagrange::CurveInterpolant chi4(p1, p3, ny);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::z,   &chi1, BC::Neumann);
            m.boundary(Field::z,   &chi2, BC::Dirichlet);
            m.boundary(Field::z,   &chi3, BC::Dirichlet);
            m.boundary(Field::z,   &chi4, BC::Neumann);

            m.boundary(Field::v1,  &chi1, BC::Neumann);
            m.boundary(Field::n11, &chi2, BC::Dirichlet, 1);
            m.boundary(Field::n12, &chi3, BC::Dirichlet);
            m.boundary(Field::v1,  &chi4, BC::Dirichlet);

            m.boundary(Field::v2,  &chi1, BC::Dirichlet);
            m.boundary(Field::n12, &chi2, BC::Dirichlet);
            m.boundary(Field::n22, &chi3, BC::Dirichlet, 1);
            m.boundary(Field::v2,  &chi4, BC::Neumann);

            m.planeStrain();
            m.load(1);
            m.linear();

            m.output("plot/Data/Membrane/z",   Field::z);
            m.output("plot/Data/Membrane/v1",  Field::v1);
            m.output("plot/Data/Membrane/v2",  Field::v2);
            m.principalStresses("plot/Data/Membrane/kartesianStresses", "plot/Data/Membrane/principalStresses");
            m.principalStrains("plot/Data/Membrane/kartesianDeformations", "plot/Data/Membrane/kartesianStrains", "plot/Data/Membrane/principalStrains");
            break;
        }
    }
}