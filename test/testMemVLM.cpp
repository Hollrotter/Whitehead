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

            arma::vec x1A = Chebyshev::gaussLobatto(30)/2;
            arma::vec x2A = (1+Chebyshev::gaussLobatto(40))/2;

            arma::vec yA = x2A;
            arma::mat xA = x1A*arma::ones(1, x2A.size());

            VLM vlm(xA, yA);
            vlm.pitch(3);
            vlm.dynamicPressure(0.1);
            vlm.symmetry(Symmetry::y);

            MemVLM memvlm(&membrane, &vlm);

            membrane.output(Field::z, "plot/Data/MemVLM/z");
            break;
        }
    }
}