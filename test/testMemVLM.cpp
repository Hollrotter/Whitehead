#include "MemVLM.hpp"

int main()
{
    switch (1)
    {
        case 0: // Rectangular Wing
        {
            Point p1(-0.5, 0);
            Point p2( 0.5, 0);
            Point p3( 0.5, 1);
            Point p4(-0.5, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, 10);
            Lagrange::CurveInterpolant chi2(p2, p3, 10);
            Lagrange::CurveInterpolant chi3(p3, p4, 10);
            Lagrange::CurveInterpolant chi4(p4, p1, 10);

            Membrane membrane({&chi1, &chi2, &chi3, &chi4});

            membrane.boundary(Field::z,   &chi1, BC::Neumann);
            membrane.boundary(Field::z,   &chi2, BC::Dirichlet);
            membrane.boundary(Field::z,   &chi3, BC::Dirichlet);
            membrane.boundary(Field::z,   &chi4, BC::Dirichlet);

            membrane.boundary(Field::v2,  &chi1, BC::Dirichlet);
            membrane.boundary(Field::n12, &chi2, BC::Dirichlet);
            membrane.boundary(Field::n22, &chi3, BC::Dirichlet, 0.15);
            membrane.boundary(Field::v2,  &chi4, BC::Neumann);

            membrane.boundary(Field::v1,  &chi1, BC::Neumann);
            membrane.boundary(Field::n11, &chi2, BC::Dirichlet, 0.15);
            membrane.boundary(Field::n12, &chi3, BC::Dirichlet);
            membrane.boundary(Field::v1,  &chi4, BC::Dirichlet);

            membrane.planeStrain();

            arma::vec x1A = Chebyshev::gaussLobatto(30)/2;
            arma::vec x2A = (1+Chebyshev::gaussLobatto(40))/2;

            arma::vec yA = x2A;
            arma::mat xA = x1A*arma::ones(1, x2A.size());

            VLM vlm(xA, yA);
            vlm.pitch(3);
            vlm.dynamicPressure(0.1);
            vlm.symmetry(Symmetry::y);

            MemVLM memvlm(&membrane, &vlm);
            memvlm(Coupling::monolithic);

            memvlm.linear();

            membrane.output(Field::z,   "plot/Data/MemVLM/z");
            membrane.output(Field::v1,  "plot/Data/MemVLM/v1");
            membrane.output(Field::v2,  "plot/Data/MemVLM/v2");
            membrane.output(Field::n11, "plot/Data/MemVLM/n11");
            membrane.output(Field::n12, "plot/Data/MemVLM/n12");
            membrane.output(Field::n22, "plot/Data/MemVLM/n22");

            vlm.output("plot/Data/MemVLM/p");

            std::cout << "cL = " << vlm.get_lift().t()   / 0.1   << '\n';
            std::cout << "cM = " << vlm.get_moment().t() / 0.1 << '\n';

            break;
        }
        case 1: // Elliptic Wing
        {
            double b = 10;
            double c = 1.39;

            size_t nxS = 10;
            size_t nyS = 20;

            size_t nxA = 40;
            size_t nyA = 60;

            arma::mat data;
            data.load(arma::csv_name("clarky-il.csv", arma::csv_opts::no_header));

            arma::vec x1A = c/2*(1+Chebyshev::gaussLobatto(nxA));
            arma::vec x2A = 0.98*b/4*(1+Chebyshev::gaussLobatto(nyA));

            arma::vec yA = x2A;
            arma::mat xA = (x1A - c/4)*sqrt(1-pow(x2A/(b/2), 2)).t();

            VLM vlm(xA, yA);
            vlm.pitch(5);
            vlm.camber(Camber(Splinefit(data.col(0)/1e3, data.col(1)/1e3, 10)));

            vlm.dynamicPressure(10);
            vlm.symmetry(Symmetry::y);

            arma::vec x1S = c/2*(1+Chebyshev::gaussLobatto(nxS));
            arma::vec x2S = 0.98*b/4*(1+Chebyshev::gaussLobatto(nyS));

            arma::mat yS = arma::ones(x1S.size(), 1)*x2S.t();
            arma::mat xS = (x1S - c/4)*sqrt(1-pow(x2S/(b/2), 2)).t();

            Lagrange::CurveInterpolant chi1(xS.col(0),         yS.col(0));
            Lagrange::CurveInterpolant chi2(xS.row(nxS-1).t(), yS.row(nxS-1).t());
            Lagrange::CurveInterpolant chi3(xS.col(nyS-1),     yS.col(nyS-1));
            Lagrange::CurveInterpolant chi4(xS.row(0).t(),     yS.row(0).t());

            Membrane membrane({&chi1, &chi2, &chi3, &chi4});
            membrane.boundary(Field::z,   Direction::N, BC::Dirichlet);
            membrane.boundary(Field::z,   Direction::S, BC::Neumann);
            membrane.boundary(Field::z,   Direction::W, BC::Dirichlet);
            membrane.boundary(Field::z,   Direction::E, BC::Dirichlet);

            membrane.boundary(Field::n12, Direction::N, BC::Dirichlet);
            membrane.boundary(Field::v1,  Direction::S, BC::Dirichlet);
            membrane.boundary(Field::v1,  Direction::W, BC::Dirichlet);
            membrane.boundary(Field::n11, Direction::E, BC::Dirichlet, 25);

            membrane.boundary(Field::n22, Direction::N, BC::Dirichlet, 25);
            membrane.boundary(Field::v2,  Direction::S, BC::Dirichlet);
            membrane.boundary(Field::n12, Direction::W, BC::Dirichlet);
            membrane.boundary(Field::n12, Direction::E, BC::Dirichlet);

            membrane.planeStrain();

            MemVLM memvlm(&membrane, &vlm);
            memvlm(Coupling::monolithic);
            memvlm.linear();

            membrane.output(Field::z,   "plot/Data/MemVLM/z");
            membrane.output(Field::v1,  "plot/Data/MemVLM/v1");
            membrane.output(Field::v2,  "plot/Data/MemVLM/v2");
            membrane.output(Field::n11, "plot/Data/MemVLM/n11");
            membrane.output(Field::n12, "plot/Data/MemVLM/n12");
            membrane.output(Field::n22, "plot/Data/MemVLM/n22");
            vlm.output("plot/Data/MemVLM/p");
            break;
        }
    }
}