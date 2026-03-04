#include "DVM.hpp"

int main()
{
    switch (1)
    {
        case 0: // Verification of linear DVM with some analytical solutions
        {
            int N = 50;
            double c = 1.0;
            double alpha = 0.1;
            double h = 0.1;

            arma::vec x(N+1, 1, arma::fill::zeros);
            arma::vec  xg(N, 1, arma::fill::zeros);
            
            for (int n = 0; n < N+1; n++)
                x(n)  = c/2*(1 - cos(arma::datum::pi*n/(N+1)));
            for (int n = 0; n < N; n++)
                xg(n) = x(n) + 0.25*(x(n+1)-x(n));

            // 1. Birnbaum-Ackermannsche Normalverteilung
            arma::vec dcp1 = 4*alpha*sqrt((1-xg/c)/(xg/c));

            // 2. Birnbaum-Ackermannsche Normalverteilung
            arma::vec dcp2 = 32*h/c*sqrt((xg/c)%(1-xg/c));

            DVM dvm1(c, N);
            dvm1.pitch(180/arma::datum::pi * alpha);
            dvm1.dvm();
            dvm1.output("plot/Data/DVM/dvm1");

            DVM dvm2(c, N, [&](double xc){ return 4*h/c*(1 - 2*xc); });
            dvm2.dvm();
            dvm2.output("plot/Data/DVM/dvm2");

            std::ofstream file("plot/Data/DVM/dcp");
            for (int n = 0; n < N; n++)
                file << xg(n) << ' ' << dcp1(n) << ' ' << dcp2(n) << '\n';
            file.close();
        
            break;
        }
        case 1: // A test for the nonlinear DVM
        {
            double c = 2.0;
            double h = 1;
            DVM dvm(c, 100, [&](double xc){ return 4*h*xc*(1 -   xc); },
                            [&](double xc){ return 4*h/c *(1 - 2*xc); });
            dvm.pitch(5);
            dvm.dynamicPressure(1.225/2*15*15);
            dvm.aerodynamics(Analysis::nonlinear);
            dvm.dvm();
            dvm.output("plot/Data/DVM/dvm");

            std::cout << "cL = " << dvm.get_lift().t()   << '\n';
            std::cout << "cM = " << dvm.get_moment().t() << '\n';
        }
    }
}