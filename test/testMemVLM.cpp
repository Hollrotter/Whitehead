#include "MemVLM.hpp"

int main()
{
    switch (0)
    {
        case 0:
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
    }
}