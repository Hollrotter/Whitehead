#include "MemVLM.hpp"

int main()
{
    size_t n1S = 6;
    size_t n2S = 6;
    size_t n1A = 6;
    size_t n2A = 6;

    std::ofstream file("plot/Data/SquareWing");
    for (size_t i = 0; i < 11; i++, n1S+=2, n2S+=2, n1A+=12, n2A+=12)
    {
        auto t1 = std::chrono::high_resolution_clock::now();

        double q_dyn = 1;

        Point p1(-1, 0);
        Point p2( 1, 0);
        Point p3( 1, 2);
        Point p4(-1, 2);
        
        std::cout << "n1S = " << n1S << " n2S = " << n2S << std::endl;
        std::cout << "n1A = " << n1A << " n2A = " << n2A << std::endl;

        Lagrange::CurveInterpolant chi1(p1, p2, n1S);
        Lagrange::CurveInterpolant chi2(p2, p3, n2S);
        Lagrange::CurveInterpolant chi3(p3, p4, n1S);
        Lagrange::CurveInterpolant chi4(p4, p1, n2S);

        Membrane membrane({&chi1, &chi2, &chi3, &chi4});

        membrane.boundary(Field::z,   &chi1, BC::Neumann);
        membrane.boundary(Field::z,   &chi2, BC::Dirichlet);
        membrane.boundary(Field::z,   &chi3, BC::Dirichlet);
        membrane.boundary(Field::z,   &chi4, BC::Dirichlet);

        membrane.boundary(Field::v1,  &chi1, BC::Neumann);
        membrane.boundary(Field::n11, &chi2, BC::Dirichlet, 6);
        membrane.boundary(Field::n12, &chi3, BC::Dirichlet);
        membrane.boundary(Field::v1,  &chi4, BC::Dirichlet);

        membrane.boundary(Field::v2,  &chi1, BC::Dirichlet);
        membrane.boundary(Field::n12, &chi2, BC::Dirichlet);
        membrane.boundary(Field::n22, &chi3, BC::Dirichlet, 6);
        membrane.boundary(Field::v2,  &chi4, BC::Neumann);

        membrane.planeStrain();

        arma::vec x1A = Chebyshev::gaussLobatto(n1A);
        arma::vec x2A = Chebyshev::gaussLobatto(n2A)+1;

        arma::vec yA = x2A;
        arma::mat xA = x1A*arma::ones(1, x2A.size());

        VLM vlm(xA, yA);
        vlm.pitch(3);
        vlm.dynamicPressure(q_dyn);
        vlm.symmetry(Symmetry::y);

        MemVLM memvlm(&membrane, &vlm);
        memvlm(Coupling::monolithic);

        memvlm.linear();
        
        double cL = vlm.get_lift()(0)   / 4 / q_dyn;
        double cM = vlm.get_moment()(0) / 8 / q_dyn;

        std::cout << "cL = " << cL << '\n';
        std::cout << "cM = " << cM << '\n';

        auto t2 = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> dt = t2 - t1;
        file << n1S << ' ' << n2S << ' ' << n1A << ' ' << n2A << ' ' << cL << ' ' << cM << ' ' << dt.count() << '\n';
    }
    file.close();
}