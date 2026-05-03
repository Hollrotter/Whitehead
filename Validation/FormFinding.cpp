#include "Structure.hpp"

int main()
{
    switch (6)
    {
        case 0: // An inclined plane
        {
            size_t n = 20;

            Point p1(-5,-5);
            Point p2( 5,-5);
            Point p3( 5, 5);
            Point p4(-5, 5);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::v1,  Direction::S, BC::Neumann);
            m.boundary(Field::n11, Direction::E, BC::Dirichlet, 1);
            m.boundary(Field::n12, Direction::N, BC::Dirichlet);
            m.boundary(Field::v1,  Direction::W, BC::Dirichlet);

            m.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            m.boundary(Field::n12, Direction::E, BC::Dirichlet);
            m.boundary(Field::n22, Direction::N, BC::Dirichlet, 1);
            m.boundary(Field::v2,  Direction::W, BC::Neumann);

            m.planeStrain();

            arma::vec north = 5*(Chebyshev::gaussLobatto(n) + 1);
            arma::vec south = 5*(Chebyshev::gaussLobatto(n) - 1);
            arma::vec west  = 5*(Chebyshev::gaussLobatto(n) - 1);
            arma::vec east  = 5*(Chebyshev::gaussLobatto(n) + 1);

            m.boundary(Field::z, Direction::S, BC::Dirichlet, south);
            m.boundary(Field::z, Direction::E, BC::Dirichlet, east);
            m.boundary(Field::z, Direction::N, BC::Dirichlet, north);
            m.boundary(Field::z, Direction::W, BC::Dirichlet, west);

            m.iterations(100);
            m.linear();

            m.boundary(Field::n11, Direction::E, BC::None);
            m.boundary(Field::n12, Direction::N, BC::None);
            m.boundary(Field::n12, Direction::E, BC::None);
            m.boundary(Field::n22, Direction::N, BC::None);

            m.boundary(Field::v1, Direction::S, BC::Dirichlet);
            m.boundary(Field::v1, Direction::E, BC::Dirichlet);
            m.boundary(Field::v1, Direction::N, BC::Dirichlet);
            m.boundary(Field::v1, Direction::W, BC::Dirichlet);

            m.boundary(Field::v2, Direction::S, BC::Dirichlet);
            m.boundary(Field::v2, Direction::E, BC::Dirichlet);
            m.boundary(Field::v2, Direction::N, BC::Dirichlet);
            m.boundary(Field::v2, Direction::W, BC::Dirichlet);

            m.nonlinear();

            m.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 1: // A tent-like structure
        {
            size_t n = 10;
            size_t m = 25;

            Point p1(-2, 0);
            Point p2( 0, 0);
            Point p3( 0, 5);
            Point p4(-2, 5);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, m);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, m);

            Membrane membrane({&chi1, &chi2, &chi3, &chi4});

            membrane.boundary(Field::v1,  Direction::S, BC::Neumann);
            membrane.boundary(Field::v1,  Direction::E, BC::Dirichlet);
            membrane.boundary(Field::n12, Direction::N, BC::Dirichlet);
            membrane.boundary(Field::n11, Direction::W, BC::Dirichlet, 1);

            membrane.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            membrane.boundary(Field::v2,  Direction::E, BC::Neumann);
            membrane.boundary(Field::n22, Direction::N, BC::Dirichlet, 1);
            membrane.boundary(Field::n12, Direction::W, BC::Dirichlet);

            membrane.planeStrain();

            arma::vec north = 0.5*(Chebyshev::gaussLobatto(n) + 1);

            membrane.boundary(Field::z, Direction::S, BC::Neumann);
            membrane.boundary(Field::z, Direction::E, BC::Neumann);
            membrane.boundary(Field::z, Direction::N, BC::Dirichlet, north);
            membrane.boundary(Field::z, Direction::W, BC::Dirichlet);

            membrane.iterations(50);
            membrane.linear();

            membrane.boundary(Field::n12, Direction::N, BC::None);
            membrane.boundary(Field::n11, Direction::W, BC::None);
            membrane.boundary(Field::n22, Direction::N, BC::None);
            membrane.boundary(Field::n12, Direction::W, BC::None);

            membrane.boundary(Field::v1,  Direction::N, BC::Dirichlet);
            membrane.boundary(Field::v1,  Direction::W, BC::Dirichlet);
            membrane.boundary(Field::v1,  Direction::S, BC::Dirichlet);
            membrane.boundary(Field::v1,  Direction::E, BC::Dirichlet);

            membrane.boundary(Field::v2,  Direction::N, BC::Dirichlet);
            membrane.boundary(Field::v2,  Direction::W, BC::Dirichlet);
            membrane.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            membrane.boundary(Field::v2,  Direction::E, BC::Dirichlet);

            membrane.nonlinear();

            membrane.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 2: // An analytical boundary
        {
            /**
             * This one is not working properly. The reason is, that z is zero at most
             * of the boundary. Only at two corners, there is a concentrated increase of z.
             * This does not make it suitable for Spectral Methods which rely on global
             * distributions rather than concentrated disturbances. We split the geometry
             * up and solve an instance of the Structure class. Even in this case, the
             * nonlinear analysis fails. We are only able to achieve good results with the
             * lienar analysis.
             */

            size_t n = 10;

            Point p1(-4,-4);
            Point p2(-2.5,-4);
            Point p3(-2.5,-2.5);
            Point p4(-4,-2.5);

            Point p5( 4, 4);
            Point p6( 2.5, 4);
            Point p7( 2.5, 2.5);
            Point p8( 4, 2.5);

            Point p9 ( 4,-4);
            Point p10( 2.5,-4);
            Point p11( 2.5,-2.5);
            Point p12( 4,-2.5);

            Point p13(-4, 4);
            Point p14(-2.5, 4);
            Point p15(-2.5, 2.5);
            Point p16(-4, 2.5);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            Lagrange::CurveInterpolant chi5(p5, p6, n);
            Lagrange::CurveInterpolant chi6(p6, p7, n);
            Lagrange::CurveInterpolant chi7(p7, p8, n);
            Lagrange::CurveInterpolant chi8(p8, p5, n);

            Lagrange::CurveInterpolant chi9 (p9,  p10, n);
            Lagrange::CurveInterpolant chi10(p10, p11, n);
            Lagrange::CurveInterpolant chi11(p11, p12, n);
            Lagrange::CurveInterpolant chi12(p12, p9,  n);

            Lagrange::CurveInterpolant chi13(p13, p14, n);
            Lagrange::CurveInterpolant chi14(p14, p15, n);
            Lagrange::CurveInterpolant chi15(p15, p16, n);
            Lagrange::CurveInterpolant chi16(p16, p13, n);

            Lagrange::CurveInterpolant chi17(p3,  p11, n);
            Lagrange::CurveInterpolant chi18(p3,  p15, n);
            Lagrange::CurveInterpolant chi19(p4,  p16, n);
            Lagrange::CurveInterpolant chi20(p2,  p10, n);

            Lagrange::CurveInterpolant chi21(p7,  p11, n);
            Lagrange::CurveInterpolant chi22(p7,  p15, n);
            Lagrange::CurveInterpolant chi23(p6,  p14, n);
            Lagrange::CurveInterpolant chi24(p8,  p12, n);

            Membrane m1({&chi1,  &chi2,  &chi3,  &chi4});
            Membrane m2({&chi7,  &chi8,  &chi5,  &chi6});
            Membrane m3({&chi9,  &chi12, &chi11, &chi10});
            Membrane m4({&chi15, &chi14, &chi13, &chi16});
            Membrane m5({&chi17, &chi21, &chi22, &chi18});
            Membrane m6({&chi20, &chi10, &chi17, &chi2});
            Membrane m7({&chi11, &chi24, &chi7,  &chi21});
            Membrane m8({&chi22, &chi6,  &chi23, &chi14});
            Membrane m9({&chi3,  &chi18, &chi15, &chi19});

            Structure s({&m1, &m2, &m3, &m4, &m5, &m6, &m7, &m8, &m9});
            s.setIterations(100);

            s.boundary(Field::v1,  &chi1,  BC::Neumann);
            s.boundary(Field::v1,  &chi20, BC::Neumann);
            s.boundary(Field::v1,  &chi9,  BC::Neumann);
            s.boundary(Field::n12, &chi13, BC::Dirichlet);
            s.boundary(Field::n12, &chi23, BC::Dirichlet);
            s.boundary(Field::n12, &chi5,  BC::Dirichlet);
            s.boundary(Field::n11, &chi4,  BC::Dirichlet, 0.1);
            s.boundary(Field::n11, &chi19, BC::Dirichlet, 0.1);
            s.boundary(Field::n11, &chi16, BC::Dirichlet, 0.1);
            s.boundary(Field::v1,  &chi12, BC::Dirichlet);
            s.boundary(Field::v1,  &chi24, BC::Dirichlet);
            s.boundary(Field::v1,  &chi8,  BC::Dirichlet);

            s.boundary(Field::v2,  &chi1,  BC::Dirichlet);
            s.boundary(Field::v2,  &chi20, BC::Dirichlet);
            s.boundary(Field::v2,  &chi9,  BC::Dirichlet);
            s.boundary(Field::n22, &chi13, BC::Dirichlet, 0.1);
            s.boundary(Field::n22, &chi23, BC::Dirichlet, 0.1);
            s.boundary(Field::n22, &chi5,  BC::Dirichlet, 0.1);
            s.boundary(Field::n12, &chi4,  BC::Dirichlet);
            s.boundary(Field::n12, &chi19, BC::Dirichlet);
            s.boundary(Field::n12, &chi16, BC::Dirichlet);
            s.boundary(Field::v2,  &chi12, BC::Neumann);
            s.boundary(Field::v2,  &chi24, BC::Neumann);
            s.boundary(Field::v2,  &chi8,  BC::Neumann);

            s.planeStrain();

            arma::vec x1 = (p2.X()-p1.X())*(1+Chebyshev::gaussLobatto(n))/2+p1.X();
            arma::vec north1 = 5*((exp(pow(x1-4, 2)+1) - exp(-pow(x1-4, 2)-1))/(exp(pow(x1-4, 2)+1) + exp(-pow(x1-4, 2)-1)));
            arma::vec south1 = 5*((exp(pow(x1+4, 2)+1) - exp(-pow(x1+4, 2)-1))/(exp(pow(x1+4, 2)+1) + exp(-pow(x1+4, 2)-1)));

            arma::vec x2 = (p10.X()-p2.X())*(1+Chebyshev::gaussLobatto(n))/2+p2.X();
            arma::vec north2 = 5*((exp(pow(x2-4, 2)+1) - exp(-pow(x2-4, 2)-1))/(exp(pow(x2-4, 2)+1) + exp(-pow(x2-4, 2)-1)));
            arma::vec south2 = 5*((exp(pow(x2+4, 2)+1) - exp(-pow(x2+4, 2)-1))/(exp(pow(x2+4, 2)+1) + exp(-pow(x2+4, 2)-1)));

            arma::vec x3 = (p9.X()-p10.X())*(1+Chebyshev::gaussLobatto(n))/2+p10.X();
            arma::vec north3 = 5*((exp(pow(x3-4, 2)+1) - exp(-pow(x3-4, 2)-1))/(exp(pow(x3-4, 2)+1) + exp(-pow(x3-4, 2)-1)));
            arma::vec south3 = 5*((exp(pow(x3+4, 2)+1) - exp(-pow(x3+4, 2)-1))/(exp(pow(x3+4, 2)+1) + exp(-pow(x3+4, 2)-1)));

            arma::vec y1 = (p4.Y()-p1.Y())*(1+Chebyshev::gaussLobatto(n))/2+p1.Y();
            arma::vec west1 = 5*((exp(pow(-4-y1, 2)+1) - exp(-pow(-4-y1, 2)-1))/(exp(pow(-4-y1, 2)+1) + exp(-pow(-4-y1, 2)-1)));
            arma::vec east1 = 5*((exp(pow(4-y1, 2)+1) - exp(-pow(4-y1, 2)-1))/(exp(pow(4-y1, 2)+1) + exp(-pow(4-y1, 2)-1)));

            arma::vec y2 = (p16.Y()-p4.Y())*(1+Chebyshev::gaussLobatto(n))/2+p4.Y();
            arma::vec west2 = 5*((exp(pow(-4-y2, 2)+1) - exp(-pow(-4-y2, 2)-1))/(exp(pow(-4-y2, 2)+1) + exp(-pow(-4-y2, 2)-1)));
            arma::vec east2 = 5*((exp(pow(4-y2, 2)+1) - exp(-pow(4-y2, 2)-1))/(exp(pow(4-y2, 2)+1) + exp(-pow(4-y2, 2)-1)));

            arma::vec y3 = (p13.Y()-p16.Y())*(1+Chebyshev::gaussLobatto(n))/2+p16.Y();
            arma::vec west3 = 5*((exp(pow(-4-y3, 2)+1) - exp(-pow(-4-y3, 2)-1))/(exp(pow(-4-y3, 2)+1) + exp(-pow(-4-y3, 2)-1)));
            arma::vec east3 = 5*((exp(pow(4-y3, 2)+1) - exp(-pow(4-y3, 2)-1))/(exp(pow(4-y3, 2)+1) + exp(-pow(4-y3, 2)-1)));

            s.boundary(Field::z, &chi1,  BC::Dirichlet, south1);
            s.boundary(Field::z, &chi20, BC::Dirichlet, south2);
            s.boundary(Field::z, &chi9,  BC::Dirichlet, south3);
            s.boundary(Field::z, &chi12, BC::Dirichlet, east1);
            s.boundary(Field::z, &chi24, BC::Dirichlet, east2);
            s.boundary(Field::z, &chi8,  BC::Dirichlet, east3);
            s.boundary(Field::z, &chi13, BC::Dirichlet, north1);
            s.boundary(Field::z, &chi23, BC::Dirichlet, north2);
            s.boundary(Field::z, &chi5,  BC::Dirichlet, north3);
            s.boundary(Field::z, &chi4,  BC::Dirichlet, west1);
            s.boundary(Field::z, &chi19, BC::Dirichlet, west2);
            s.boundary(Field::z, &chi16, BC::Dirichlet, west3);

            s.linear();

            s.boundary(Field::n12, &chi13, BC::None);
            s.boundary(Field::n12, &chi23, BC::None);
            s.boundary(Field::n12, &chi5,  BC::None);
            s.boundary(Field::n11, &chi4,  BC::None);
            s.boundary(Field::n11, &chi19, BC::None);
            s.boundary(Field::n11, &chi16, BC::None);
            s.boundary(Field::n22, &chi13, BC::None);
            s.boundary(Field::n22, &chi23, BC::None);
            s.boundary(Field::n22, &chi5,  BC::None);
            s.boundary(Field::n12, &chi4,  BC::None);
            s.boundary(Field::n12, &chi19, BC::None);
            s.boundary(Field::n12, &chi16, BC::None);

            s.boundary(Field::v1, &chi1,  BC::Dirichlet);
            s.boundary(Field::v1, &chi20, BC::Dirichlet);
            s.boundary(Field::v1, &chi9,  BC::Dirichlet);
            s.boundary(Field::v1, &chi13, BC::Dirichlet);
            s.boundary(Field::v1, &chi23, BC::Dirichlet);
            s.boundary(Field::v1, &chi5,  BC::Dirichlet);
            s.boundary(Field::v1, &chi4,  BC::Dirichlet);
            s.boundary(Field::v1, &chi19, BC::Dirichlet);
            s.boundary(Field::v1, &chi16, BC::Dirichlet);
            s.boundary(Field::v1, &chi12, BC::Dirichlet);
            s.boundary(Field::v1, &chi24, BC::Dirichlet);
            s.boundary(Field::v1, &chi8,  BC::Dirichlet);

            s.boundary(Field::v2, &chi1,  BC::Dirichlet);
            s.boundary(Field::v2, &chi20, BC::Dirichlet);
            s.boundary(Field::v2, &chi9,  BC::Dirichlet);
            s.boundary(Field::v2, &chi13, BC::Dirichlet);
            s.boundary(Field::v2, &chi23, BC::Dirichlet);
            s.boundary(Field::v2, &chi5,  BC::Dirichlet);
            s.boundary(Field::v2, &chi4,  BC::Dirichlet);
            s.boundary(Field::v2, &chi19, BC::Dirichlet);
            s.boundary(Field::v2, &chi16, BC::Dirichlet);
            s.boundary(Field::v2, &chi12, BC::Dirichlet);
            s.boundary(Field::v2, &chi24, BC::Dirichlet);
            s.boundary(Field::v2, &chi8,  BC::Dirichlet);

            // s.nonlinear();

            s.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 3: // Scherk
        {
            size_t n = 20;

            Point p1(-0.5,-0.5);
            Point p2( 0.5,-0.5);
            Point p3( 0.5, 0.5);
            Point p4(-0.5, 0.5);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::v1,  Direction::S, BC::Neumann);
            m.boundary(Field::v1,  Direction::E, BC::Dirichlet);
            m.boundary(Field::n12, Direction::N, BC::Dirichlet);
            m.boundary(Field::n11, Direction::W, BC::Dirichlet, 10);

            m.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            m.boundary(Field::v2,  Direction::E, BC::Neumann);
            m.boundary(Field::n22, Direction::N, BC::Dirichlet, 10);
            m.boundary(Field::n12, Direction::W, BC::Dirichlet);

            m.planeStrain();

            arma::vec x = 0.5*Chebyshev::gaussLobatto(n);
            arma::vec y = 0.5*Chebyshev::gaussLobatto(n);

            arma::vec north = log(cos( 0.5)) - log(cos(x));
            arma::vec south = log(cos(-0.5)) - log(cos(x));
            arma::vec west  = log(cos(y)) - log(cos(-0.5));
            arma::vec east  = log(cos(y)) - log(cos( 0.5));

            m.boundary(Field::z, Direction::S, BC::Dirichlet, south);
            m.boundary(Field::z, Direction::E, BC::Dirichlet, east);
            m.boundary(Field::z, Direction::N, BC::Dirichlet, north);
            m.boundary(Field::z, Direction::W, BC::Dirichlet, west);

            m.iterations(20);
            m.linear();

            m.boundary(Field::n12, Direction::N, BC::None);
            m.boundary(Field::n11, Direction::W, BC::None);
            m.boundary(Field::n22, Direction::N, BC::None);
            m.boundary(Field::n12, Direction::W, BC::None);

            m.boundary(Field::v1, Direction::S, BC::Dirichlet);
            m.boundary(Field::v1, Direction::E, BC::Dirichlet);
            m.boundary(Field::v1, Direction::N, BC::Dirichlet);
            m.boundary(Field::v1, Direction::W, BC::Dirichlet);

            m.boundary(Field::v2, Direction::S, BC::Dirichlet);
            m.boundary(Field::v2, Direction::E, BC::Dirichlet);
            m.boundary(Field::v2, Direction::N, BC::Dirichlet);
            m.boundary(Field::v2, Direction::W, BC::Dirichlet);

            m.nonlinear();

            m.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 4: // Another analytical boundary
        {
            size_t n = 15;

            Point p1(-1,-1);
            Point p2( 1,-1);
            Point p3( 1, 1);
            Point p4(-1, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::v1,  Direction::S, BC::Neumann);
            m.boundary(Field::v1,  Direction::E, BC::Dirichlet);
            m.boundary(Field::n12, Direction::N, BC::Dirichlet);
            m.boundary(Field::n11, Direction::W, BC::Dirichlet, 10);

            m.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            m.boundary(Field::v2,  Direction::E, BC::Neumann);
            m.boundary(Field::n22, Direction::N, BC::Dirichlet, 10);
            m.boundary(Field::n12, Direction::W, BC::Dirichlet);

            m.planeStrain();

            arma::vec x = Chebyshev::gaussLobatto(n);
            arma::vec y = Chebyshev::gaussLobatto(n);

            arma::vec north = 1/(1+6*exp(-pow(x-1, 2)-1));
            arma::vec south = 1/(1+6*exp(-pow(x+1, 2)-1));
            arma::vec west  = 1/(1+6*exp(-pow(-1-y, 2)-1));
            arma::vec east  = 1/(1+6*exp(-pow(1-y, 2)-1));

            m.boundary(Field::z, Direction::S, BC::Dirichlet, south);
            m.boundary(Field::z, Direction::E, BC::Dirichlet, east);
            m.boundary(Field::z, Direction::N, BC::Dirichlet, north);
            m.boundary(Field::z, Direction::W, BC::Dirichlet, west);

            m.iterations(100);
            m.linear();

            m.boundary(Field::n12, Direction::N, BC::None);
            m.boundary(Field::n11, Direction::W, BC::None);
            m.boundary(Field::n22, Direction::N, BC::None);
            m.boundary(Field::n12, Direction::W, BC::None);

            m.boundary(Field::v1, Direction::S, BC::Dirichlet);
            m.boundary(Field::v1, Direction::E, BC::Dirichlet);
            m.boundary(Field::v1, Direction::N, BC::Dirichlet);
            m.boundary(Field::v1, Direction::W, BC::Dirichlet);

            m.boundary(Field::v2, Direction::S, BC::Dirichlet);
            m.boundary(Field::v2, Direction::E, BC::Dirichlet);
            m.boundary(Field::v2, Direction::N, BC::Dirichlet);
            m.boundary(Field::v2, Direction::W, BC::Dirichlet);

            m.nonlinear();

            m.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 5: // Hypar
        {
            size_t n = 10;

            Point p1(-25,-25);
            Point p2( 25,-25);
            Point p3( 25, 25);
            Point p4(-25, 25);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, n);

            Membrane m({&chi1, &chi2, &chi3, &chi4});

            m.boundary(Field::v1,  Direction::S, BC::Neumann);
            m.boundary(Field::v1,  Direction::E, BC::Dirichlet);
            m.boundary(Field::n12, Direction::N, BC::Dirichlet);
            m.boundary(Field::n11, Direction::W, BC::Dirichlet, 10);

            m.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            m.boundary(Field::v2,  Direction::E, BC::Neumann);
            m.boundary(Field::n22, Direction::N, BC::Dirichlet, 10);
            m.boundary(Field::n12, Direction::W, BC::Dirichlet);

            m.planeStrain();

            arma::vec x = 25*Chebyshev::gaussLobatto(n);
            arma::vec y = 25*Chebyshev::gaussLobatto(n);

            arma::vec north = 0.2*(x+25);
            arma::vec south = 0.2*(25-x);
            arma::vec west  = 0.2*(25-y);
            arma::vec east  = 0.2*(y+25);

            m.boundary(Field::z, Direction::S, BC::Dirichlet, south);
            m.boundary(Field::z, Direction::E, BC::Dirichlet, east);
            m.boundary(Field::z, Direction::N, BC::Dirichlet, north);
            m.boundary(Field::z, Direction::W, BC::Dirichlet, west);

            m.iterations(20);
            m.linear();

            m.boundary(Field::n12, Direction::N, BC::None);
            m.boundary(Field::n11, Direction::W, BC::None);
            m.boundary(Field::n22, Direction::N, BC::None);
            m.boundary(Field::n12, Direction::W, BC::None);

            m.boundary(Field::v1, Direction::S, BC::Dirichlet);
            m.boundary(Field::v1, Direction::E, BC::Dirichlet);
            m.boundary(Field::v1, Direction::N, BC::Dirichlet);
            m.boundary(Field::v1, Direction::W, BC::Dirichlet);

            m.boundary(Field::v2, Direction::S, BC::Dirichlet);
            m.boundary(Field::v2, Direction::E, BC::Dirichlet);
            m.boundary(Field::v2, Direction::N, BC::Dirichlet);
            m.boundary(Field::v2, Direction::W, BC::Dirichlet);

            m.nonlinear();

            m.output(Field::z, "plot/Data/FormFinding");

            break;
        }
        case 6: // Twin hypar
        {
            size_t n = 15;
            size_t m = 15;

            Point p1(-30,-15);
            Point p2(  0,-15);
            Point p3( 30,-15);
            Point p4( 30, 15);
            Point p5(  0, 15);
            Point p6(-30, 15);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, n);
            Lagrange::CurveInterpolant chi3(p3, p4, m);
            Lagrange::CurveInterpolant chi4(p4, p5, n);
            Lagrange::CurveInterpolant chi5(p5, p6, n);
            Lagrange::CurveInterpolant chi6(p6, p1, m);
            Lagrange::CurveInterpolant chi7(p2, p5, m);

            Membrane m1({&chi1, &chi7, &chi5, &chi6});
            Membrane m2({&chi2, &chi3, &chi4, &chi7});

            Structure s({&m1, &m2});
            s.setIterations(200);

            s.boundary(Field::v1,  &chi1, BC::Neumann);
            s.boundary(Field::v1,  &chi2, BC::Neumann);
            s.boundary(Field::n11, &chi3, BC::Dirichlet, 1);
            s.boundary(Field::n12, &chi4, BC::Dirichlet);
            s.boundary(Field::n12, &chi5, BC::Dirichlet);
            s.boundary(Field::v1,  &chi6, BC::Dirichlet);

            s.boundary(Field::v2,  &chi1, BC::Dirichlet);
            s.boundary(Field::v2,  &chi2, BC::Dirichlet);
            s.boundary(Field::n12, &chi3, BC::Dirichlet);
            s.boundary(Field::n22, &chi4, BC::Dirichlet, 1);
            s.boundary(Field::n22, &chi5, BC::Dirichlet, 1);
            s.boundary(Field::v2,  &chi6, BC::Neumann);

            s.planeStrain();

            arma::vec z1 = 2.5*(1-Chebyshev::gaussLobatto(n));
            arma::vec z2 = 2.5*(1+Chebyshev::gaussLobatto(n));
            arma::vec z3 = 2.5*(1-Chebyshev::gaussLobatto(m));
            arma::vec z4 = 2.5*(1-Chebyshev::gaussLobatto(n));
            arma::vec z5 = 2.5*(1+Chebyshev::gaussLobatto(n));
            arma::vec z6 = 2.5*(1-Chebyshev::gaussLobatto(m));

            s.boundary(Field::z,   &chi1, BC::Dirichlet, z1);
            s.boundary(Field::z,   &chi2, BC::Dirichlet, z2);
            s.boundary(Field::z,   &chi3, BC::Dirichlet, z3);
            s.boundary(Field::z,   &chi4, BC::Dirichlet, z4);
            s.boundary(Field::z,   &chi5, BC::Dirichlet, z5);
            s.boundary(Field::z,   &chi6, BC::Dirichlet, z6);

            s.linear();

            s.boundary(Field::n11, &chi3, BC::None);
            s.boundary(Field::n12, &chi4, BC::None);
            s.boundary(Field::n12, &chi5, BC::None);
            s.boundary(Field::n12, &chi3, BC::None);
            s.boundary(Field::n22, &chi4, BC::None);
            s.boundary(Field::n22, &chi5, BC::None);

            s.boundary(Field::v1,  &chi1, BC::Dirichlet);
            s.boundary(Field::v1,  &chi2, BC::Dirichlet);
            s.boundary(Field::v1,  &chi3, BC::Dirichlet);
            s.boundary(Field::v1,  &chi4, BC::Dirichlet);
            s.boundary(Field::v1,  &chi5, BC::Dirichlet);
            s.boundary(Field::v1,  &chi6, BC::Dirichlet);

            s.boundary(Field::v2,  &chi1, BC::Dirichlet);
            s.boundary(Field::v2,  &chi2, BC::Dirichlet);
            s.boundary(Field::v2,  &chi3, BC::Dirichlet);
            s.boundary(Field::v2,  &chi4, BC::Dirichlet);
            s.boundary(Field::v2,  &chi5, BC::Dirichlet);
            s.boundary(Field::v2,  &chi6, BC::Dirichlet);

            s.nonlinear();

            s.output(Field::z, "plot/Data/FormFinding");
            s.principalStrains("plot/Data/FormFinding_vx", "plot/Data/FormFinding_gx", "plot/Data/FormFinding_g1");

            break;
        }
    }
}