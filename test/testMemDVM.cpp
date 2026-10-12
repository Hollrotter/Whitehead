#include "MemDVM.hpp"

int main()
{
    double alpha = 3;
    double c = 1.0;
    double v_inf = 15;
    double rho_inf = 1.225;
    double q_inf = rho_inf/2*pow(v_inf, 2);
    double CT = 3;
    double T = CT*q_inf*c;

    DVM dvm(c, 30);
    dvm.dynamicPressure(q_inf);
    dvm.pitch(alpha);

    String string(c, 20, T);
    string.boundary(Location::Front, BC::Dirichlet);
    string.boundary(Location::Back,  BC::Dirichlet);

    MemDVM memdvm(&dvm, &string);
    memdvm(Coupling::monolithic);
    memdvm.linear();
    string.output("plot/Data/MemDVM/z");
    dvm.output("plot/Data/MemDVM/p");

    std::cout << "cL = " << dvm.get_lift().t()   << '\n';
    std::cout << "cM = " << dvm.get_moment().t() << '\n';
}