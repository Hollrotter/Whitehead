#include "Airfoil.hpp"

int main()
{
    switch (1)
    {
        case 0: // Flat plate
        {
            double c = 2;
            size_t N = 20;
            double alpha = 0.1;

            // 1. Birnbaum-Ackermannsche Normalverteilung
            arma::vec x = c/2*(1 + Chebyshev::gauss(N));
            arma::vec dcp = 4*alpha*sqrt((1-x/c)/(x/c));

            Airfoil airfoil(c, N);
            airfoil.pitch(180/arma::datum::pi * alpha);
            airfoil.linear();
            airfoil.output("plot/Data/Airfoil/flat");

            std::ofstream file("plot/Data/Airfoil/dcp");
            for (int n = 0; n < N; n++)
                file << x(n) << ' ' << dcp(n) << '\n';
            file.close();

            std::cout << "cL = " << airfoil.get_lift().t()   << std::endl;
            std::cout << "cM = " << airfoil.get_moment().t() << std::endl;
            break;
        }
        case 1: // Nonlinear
        {
            double c = 2;
            size_t n = 10;
            double h = 0.1;
            arma::vec s = Chebyshev::gauss(n);
            arma::vec x = c/2*(1+s);
            arma::vec z = h*(1-pow(s, 2));
            Lagrange::CurveInterpolant chi(x, z, s);
            Airfoil airfoil(&chi);
            airfoil.pitch(2);
            airfoil.nonlinear();
            airfoil.output("plot/Data/Airfoil/nonlinear");
            break;
        }
    }
}