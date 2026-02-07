#include "Metric.hpp"
#include "Chebyshev.hpp"

int main()
{
    double c = 1.39;
    double b = 10;
    double eps = 0.6;

    size_t n = 15;
    size_t m = 15;

    arma::mat Ix = arma::ones(n, 1);
    arma::mat Iy = arma::ones(1, m);

    arma::vec x1 = Chebyshev::gaussLobatto(n);
    arma::vec x2 = Chebyshev::gaussLobatto(m);
    arma::mat D1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
    arma::mat D2 = Chebyshev::derivativeMatrix(x2, Derivative::first);
    arma::mat D11 = Chebyshev::derivativeMatrix(x1, Derivative::second);
    arma::mat D22 = Chebyshev::derivativeMatrix(x2, Derivative::second);
    arma::vec x1_s = c*(1+x1)/2;
    arma::vec x2_s = eps*b/2*(1+x2)/2;
    arma::mat y = arma::ones(x1_s.size(), 1) * x2_s.t();
    arma::mat x = (x1_s-c/4) * sqrt(1-pow(x2_s/(b/2), 2)).t();
    auto [g_c, gc] = Metric(x, y, D1, D2);
    
    arma::mat g_11 = g_c.slice(0);
    arma::mat g_12 = g_c.slice(1);
    arma::mat g_22 = g_c.slice(2);
    arma::mat g11 = gc.slice(0);
    arma::mat g12 = gc.slice(1);
    arma::mat g22 = gc.slice(2);

    arma::mat g_11_a = -c/4*c/4*Ix * (eps*eps*pow(x2 + 1, 2) - 4).t();
    arma::mat g_12_a = -(c/4*c/4*eps*eps*(2*x1 + 1) * (2*x2 + 2).t())/4;
    arma::mat g_22_a = (b/2*b/2*eps*eps)/4 - (c/4*c/4*pow(eps, 4)*pow(2*x1 + 1, 2)*pow(x2 + 1, 2).t())/(4*(eps*eps*pow(Ix*x2.t() + 1, 2) - 4));
    arma::mat g11_a = -(4*((b/2*b/2*eps*eps)/4 - (c/4*c/4*pow(eps, 4)*pow(2*x1 + 1, 2)*pow(x2 + 1, 2).t())/(4*(eps*eps*pow(Ix*x2.t() + 1, 2) - 4))))/(b/2*b/2*c/4*c/4*eps*eps*(eps*eps*Ix*(x2%x2).t() + 2*eps*eps*Ix*x2.t() + eps*eps - 4));
    arma::mat g12_a = -((2*x1 + 1) * (((2*x2 + 2))/(b/2*b/2*(eps*eps*x2%x2 + 2*eps*eps*x2 + eps*eps - 4))).t());
    arma::mat g22_a = 4/pow(b/2*eps, 2)*arma::ones(n, m);

    std::cout << "g_11\n" << g_11 - g_11_a << "\n\n";
    std::cout << "g_12\n" << g_12 - g_12_a << "\n\n";
    std::cout << "g_22\n" << g_22 - g_22_a << "\n\n";
    std::cout << "g11\n"  << g11  - g11_a  << "\n\n";
    std::cout << "g12\n"  << g12  - g12_a  << "\n\n";
    std::cout << "g22\n"  << g22  - g22_a  << "\n\n";

    arma::mat dg11d1 = D1*g_11;
    arma::mat dg12d1 = D1*g_12;
    arma::mat dg22d1 = D1*g_22;
    arma::mat dg11d2 = g_11*D2.t();
    arma::mat dg12d2 = g_12*D2.t();
    arma::mat dg22d2 = g_22*D2.t();

    arma::mat dg11d1_a = arma::zeros(n, m);
    arma::mat dg12d1_a =-(c/4*c/4*eps*eps*(2*Ix*x2.t() + 2))/2;
    arma::mat dg22d1_a =-(c/4*c/4*pow(eps, 4)*(8*x1 + 4) * ((pow(x2 + 1, 2))/(4*(eps*eps*pow(x2 + 1, 2) - 4))).t());
    arma::mat dg11d2_a =-c/4*c/4*eps*eps*(2*Ix*x2.t() + 2);
    arma::mat dg12d2_a =-(c/4*c/4*eps*eps*(2*x1*Iy + 1))/2;
    arma::mat dg22d2_a = (2*c/4*c/4*pow(eps, 4)*pow(2*x1 + 1, 2) * (((x2 + 1))/pow(eps*eps*x2%x2 + 2*eps*eps*x2 + eps*eps - 4, 2)).t());

    std::cout << "dg11d1\n"  << dg11d1  - dg11d1_a  << "\n\n";
    std::cout << "dg12d1\n"  << dg12d1  - dg12d1_a  << "\n\n";
    std::cout << "dg22d1\n"  << dg22d1  - dg22d1_a  << "\n\n";
    std::cout << "dg11d2\n"  << dg11d2  - dg11d2_a  << "\n\n";
    std::cout << "dg12d2\n"  << dg12d2  - dg12d2_a  << "\n\n";
    std::cout << "dg22d2\n"  << dg22d2  - dg22d2_a  << "\n\n";

    arma::cube gam = Christoffel(g_c, gc, D1, D2);
    arma::mat gam111 = gam.slice(0);
    arma::mat gam112 = gam.slice(1);
    arma::mat gam122 = gam.slice(2);
    arma::mat gam211 = gam.slice(3);
    arma::mat gam212 = gam.slice(4);
    arma::mat gam222 = gam.slice(5);

    arma::mat gam111_a = arma::zeros(n, m);
    arma::mat gam112_a = Ix * ((eps*eps*(x2 + 1))/(pow(eps*x2, 2) + 2*eps*eps*x2 + eps*eps - 4)).t();
    arma::mat gam122_a = -(2*eps*eps*(2*x1 + 1)) * (1/pow(pow(eps*x2, 2) + 2*eps*eps*x2 + eps*eps - 4, 2)).t();
    arma::mat gam211_a = arma::zeros(n, m);
    arma::mat gam212_a = arma::zeros(n, m);
    arma::mat gam222_a = arma::zeros(n, m);

    std::cout << "g111\n" << gam111 - gam111_a << "\n\n";
    std::cout << "g112\n" << gam112 - gam112_a << "\n\n";
    std::cout << "g122\n" << gam122 - gam122_a << "\n\n";
    std::cout << "g211\n" << gam211 - gam211_a << "\n\n";
    std::cout << "g212\n" << gam212 - gam212_a << "\n\n";
    std::cout << "g222\n" << gam222 - gam222_a << "\n\n";
}