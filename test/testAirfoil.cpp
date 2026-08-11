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
            double h = 0.1;

            arma::vec x = c/2*(1 + Chebyshev::gauss(5*N));

            // 1. Birnbaum-Ackermannsche Normalverteilung
            arma::vec dcp1 = 4*alpha*sqrt((1-x/c)/(x/c));

            // 2. Birnbaum-Ackermannsche Normalverteilung
            arma::vec dcp2 = 32*h/c*sqrt((x/c)%(1-x/c));

            Airfoil airfoil1(c, N);
            airfoil1.pitch(180/arma::datum::pi * alpha);
            airfoil1.linear();
            airfoil1.output("plot/Data/Airfoil/flat1");

            printf("dcLdalpha_exact     = %.10f\n", arma::datum::tau);
            printf("dcLdalpha_numerical = %.10f\n", airfoil1.get_lift()/alpha);

            Airfoil airfoil2(c, N, [&](double xc){ return 4*h/c*(1 - 2*xc); });
            airfoil2.linear();
            airfoil2.output("plot/Data/Airfoil/flat2");
            
            printf("cL_exact     = %.10f\n", arma::datum::tau);
            printf("cL_numerical = %.10f\n", airfoil2.get_lift()/(2*h/c));

            std::ofstream file("plot/Data/Airfoil/dcp");
            for (size_t n = 0; n < 5*N; n++)
                file << x(n) << ' ' << dcp1(n) << ' ' << dcp2(n) << '\n';
            file.close();

            break;
        }
        case 1: // Nonlinear
        {
            double c = 1;
            size_t n = 19;
            double h = 0.5;
            arma::vec s = Chebyshev::gauss(n);
            arma::vec x = c/2*(1+s);
            arma::vec z = h*(1-pow(s, 2));
            Lagrange::CurveInterpolant chi(x, z, s);
            Airfoil airfoil(&chi);
            airfoil.pitch(10);
            airfoil.nonlinear();
            airfoil.output("plot/Data/Airfoil/nonlinear");

            std::cout << "cL = " << airfoil.get_lift()   << '\n';
            std::cout << "cM = " << airfoil.get_moment() << '\n';
            break;
        }
    }
}