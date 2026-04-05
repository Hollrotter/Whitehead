#include <ginac/ginac.h>
#include "Metric.hpp"
#include "Chebyshev.hpp"
#include "Lagrange.hpp"

using namespace GiNaC;

int main()
{
    symbol x1("x1"), x2("x2");
    const size_t samples = 1;
    size_t nx = 12;
    size_t ny = 12;
    double pi = arma::datum::pi;

    // Geometry

    // ex x = x1;
    // ex y = x2;

    double a1 = 0.2;
    double a2 = 0.1;

    ex x = x1 + a1*sin(pi*x2)*cos(pi/8*x1);
    ex y = x2 + a2*sin(pi*x1)*cos(pi/8*x2);
    
    // ex x = (x1+2)*cos(x2);
    // ex y = (x1+2)*sin(x2);

    // double b = 10;
    // double c = 1.39;

    // ex x1_a = c/2*(1+x1);
    // ex x2_a = 0.9*b/4*(1+x2);

    // ex y = x2_a;
    // ex x = (x1_a-c/4)*sqrt(1-pow(x2_a/(b/2), 2));

    ex dxdx1 = diff(x, x1);
    ex dxdx2 = diff(x, x2);
    ex dydx1 = diff(y, x1);
    ex dydx2 = diff(y, x2);
    
    ex e_11 = dxdx1*dxdx1 + dydx1*dydx1;
    ex e_12 = dxdx1*dxdx2 + dydx1*dydx2;
    ex e_22 = dxdx2*dxdx2 + dydx2*dydx2;
    
    std::cout << "e_11 = \n" << e_11 << '\n';
    std::cout << "e_12 = \n" << e_12 << '\n';
    std::cout << "e_22 = \n" << e_22 << '\n';

    ex e = e_11*e_22 - e_12*e_12;
    
    ex e11 = e_22/e;
    ex e12 =-e_12/e;
    ex e22 = e_11/e;
    
    std::cout << "e11 = \n" << e11 << '\n';
    std::cout << "e12 = \n" << e12 << '\n';
    std::cout << "e22 = \n" << e22 << '\n';

    ex de_11dx1 = diff(e_11, x1);
    ex de_12dx1 = diff(e_12, x1);
    ex de_22dx1 = diff(e_22, x1);
    ex de_11dx2 = diff(e_11, x2);
    ex de_12dx2 = diff(e_12, x2);
    ex de_22dx2 = diff(e_22, x2);
    
    ex g111 = (e11*de_11dx1 + e12*(2*de_12dx1 - de_11dx2))/2;
    ex g112 = (e11*de_11dx2 + e12*de_22dx1)/2;
    ex g122 = (e11*(2*de_12dx2 - de_22dx1) + e12*de_22dx2)/2;
    ex g211 = (e12*de_11dx1 + e22*(2*de_12dx1 - de_11dx2))/2;
    ex g212 = (e12*de_11dx2 + e22*de_22dx1)/2;
    ex g222 = (e12*(2*de_12dx2 - de_22dx1) + e22*de_22dx2)/2;

    std::cout << "g111 = \n" << g111 << '\n';
    std::cout << "g112 = \n" << g112 << '\n';
    std::cout << "g122 = \n" << g122 << '\n';
    std::cout << "g211 = \n" << g211 << '\n';
    std::cout << "g212 = \n" << g212 << '\n';
    std::cout << "g222 = \n" << g222 << '\n';

    std::ofstream file("plot/Data/MetricMMS");

    for (size_t k = 0; k < samples; k++, nx+=2, ny+=2)
    {
        arma::vec x1s = Chebyshev::gaussLobatto(nx);
        arma::vec x2s = Chebyshev::gaussLobatto(ny);

        arma::mat ys(nx, ny);
        arma::mat xs(nx, ny);
        arma::mat dxdx1_A(nx, ny);
        arma::mat dxdx2_A(nx, ny);
        arma::mat dydx1_A(nx, ny);
        arma::mat dydx2_A(nx, ny);
        arma::vec e_11_A(nx*ny);
        arma::vec e_12_A(nx*ny);
        arma::vec e_22_A(nx*ny);
        arma::vec e11_A(nx*ny);
        arma::vec e12_A(nx*ny);
        arma::vec e22_A(nx*ny);
        arma::vec g111_A(nx*ny);
        arma::vec g112_A(nx*ny);
        arma::vec g122_A(nx*ny);
        arma::vec g211_A(nx*ny);
        arma::vec g212_A(nx*ny);
        arma::vec g222_A(nx*ny);

        exmap mp;
        
        for (size_t i = 0; i < nx; i++)
        {
            mp[x1] = x1s(i);
            for (size_t j = 0; j < ny; j++)
            {
                mp[x2] = x2s(j);
                xs(i, j)       = ex_to<numeric>(evalf(x   .subs(mp))).to_double();
                ys(i, j)       = ex_to<numeric>(evalf(y   .subs(mp))).to_double();
                dxdx1_A(i, j)  = ex_to<numeric>(evalf(dxdx1.subs(mp))).to_double();
                dxdx2_A(i, j)  = ex_to<numeric>(evalf(dxdx2.subs(mp))).to_double();
                dydx1_A(i, j)  = ex_to<numeric>(evalf(dydx1.subs(mp))).to_double();
                dydx2_A(i, j)  = ex_to<numeric>(evalf(dydx2.subs(mp))).to_double();
                e_11_A(i+j*nx) = ex_to<numeric>(evalf(e_11.subs(mp))).to_double();
                e_12_A(i+j*nx) = ex_to<numeric>(evalf(e_12.subs(mp))).to_double();
                e_22_A(i+j*nx) = ex_to<numeric>(evalf(e_22.subs(mp))).to_double();
                e11_A (i+j*nx) = ex_to<numeric>(evalf(e11 .subs(mp))).to_double();
                e12_A (i+j*nx) = ex_to<numeric>(evalf(e12 .subs(mp))).to_double();
                e22_A (i+j*nx) = ex_to<numeric>(evalf(e22 .subs(mp))).to_double();
                g111_A(i+j*nx) = ex_to<numeric>(evalf(g111.subs(mp))).to_double();
                g112_A(i+j*nx) = ex_to<numeric>(evalf(g112.subs(mp))).to_double();
                g122_A(i+j*nx) = ex_to<numeric>(evalf(g122.subs(mp))).to_double();
                g211_A(i+j*nx) = ex_to<numeric>(evalf(g211.subs(mp))).to_double();
                g212_A(i+j*nx) = ex_to<numeric>(evalf(g212.subs(mp))).to_double();
                g222_A(i+j*nx) = ex_to<numeric>(evalf(g222.subs(mp))).to_double();
            }
        }

        arma::mat D1  = Chebyshev::derivativeMatrix(x1s, Derivative::first);
        arma::mat D2  = Chebyshev::derivativeMatrix(x2s, Derivative::first);
        arma::field<arma::mat> J  = Jacobian(xs,  ys,  D1, D2);
        arma::cube e_c = MetricCo(J);
        arma::cube ec  = MetricContra(e_c);
        arma::cube gam = Christoffel(e_c, ec, D1, D2);

        printf("x\n");
        std::cout << xs << '\n';
        printf("y\n");
        std::cout << ys << '\n';
        printf("J_11\n");
        std::cout << J(0, 0) << '\n';
        std::cout << dxdx1_A << '\n';
        printf("J_12\n");
        std::cout << J(0, 1) << '\n';
        std::cout << dxdx2_A << '\n';
        printf("J_21\n");
        std::cout << J(1, 0) << '\n';
        std::cout << dydx1_A << '\n';
        printf("J_22\n");
        std::cout << J(1, 1) << '\n';
        std::cout << dydx2_A << '\n';
        printf("e_11\n");
        std::cout << e_c.slice(0) << '\n';
        std::cout << arma::reshape(e_11_A, nx, ny) << '\n';
        printf("e_12\n");
        std::cout << e_c.slice(1) << '\n';
        std::cout << arma::reshape(e_12_A, nx, ny) << '\n';
        printf("e_22\n");
        std::cout << e_c.slice(2) << '\n';
        std::cout << arma::reshape(e_22_A, nx, ny) << '\n';
        printf("e11\n");
        std::cout << ec.slice(0)  << '\n';
        std::cout << arma::reshape(e11_A, nx, ny) << '\n';
        printf("e12\n"); 
        std::cout << ec.slice(1)  << '\n';
        std::cout << arma::reshape(e12_A, nx, ny) << '\n';
        printf("e22\n");
        std::cout << ec.slice(2)  << '\n';
        std::cout << arma::reshape(e22_A, nx, ny) << '\n';
        printf("g111\n");
        std::cout << gam.slice(0) << '\n';
        std::cout << arma::reshape(g111_A, nx, ny) << '\n';
        printf("g112\n");
        std::cout << gam.slice(1) << '\n';
        std::cout << arma::reshape(g112_A, nx, ny) << '\n';
        printf("g122\n");
        std::cout << gam.slice(2) << '\n';
        std::cout << arma::reshape(g122_A, nx, ny) << '\n';
        printf("g211\n");
        std::cout << gam.slice(3) << '\n';
        std::cout << arma::reshape(g211_A, nx, ny) << '\n';
        printf("g212\n");
        std::cout << gam.slice(4) << '\n';
        std::cout << arma::reshape(g212_A, nx, ny) << '\n';
        printf("g222\n");
        std::cout << gam.slice(5) << '\n';
        std::cout << arma::reshape(g222_A, nx, ny) << '\n';

        arma::vec res_e_11 = vectorise(e_c.slice(0)) - e_11_A;
        arma::vec res_e_12 = vectorise(e_c.slice(1)) - e_12_A;
        arma::vec res_e_22 = vectorise(e_c.slice(2)) - e_22_A;
        arma::vec res_e11  = vectorise(ec.slice(0))  - e11_A;
        arma::vec res_e12  = vectorise(ec.slice(1))  - e12_A;
        arma::vec res_e22  = vectorise(ec.slice(2))  - e22_A;
        arma::vec res_g111 = vectorise(gam.slice(0)) - g111_A;
        arma::vec res_g112 = vectorise(gam.slice(1)) - g112_A;
        arma::vec res_g122 = vectorise(gam.slice(2)) - g122_A;
        arma::vec res_g211 = vectorise(gam.slice(3)) - g211_A;
        arma::vec res_g212 = vectorise(gam.slice(4)) - g212_A;
        arma::vec res_g222 = vectorise(gam.slice(5)) - g222_A;
        double res_11 = sqrt(sum(res_e_11%res_e_11)/sum(e_11_A%e_11_A));
        double res_12 = sqrt(sum(res_e_12%res_e_12)/sum(e_12_A%e_12_A));
        double res_22 = sqrt(sum(res_e_22%res_e_22)/sum(e_22_A%e_22_A));
        double res11  = sqrt(sum(res_e11 %res_e11 )/sum(e11_A % e11_A));
        double res12  = sqrt(sum(res_e12 %res_e12 )/sum(e12_A % e12_A));
        double res22  = sqrt(sum(res_e22 %res_e22 )/sum(e22_A % e22_A));
        double res111 = sqrt(sum(res_g111%res_g111)/sum(g111_A%g111_A));
        double res112 = sqrt(sum(res_g112%res_g112)/sum(g112_A%g112_A));
        double res122 = sqrt(sum(res_g122%res_g122)/sum(g122_A%g122_A));
        double res211 = sqrt(sum(res_g211%res_g211)/sum(g211_A%g211_A));
        double res212 = sqrt(sum(res_g212%res_g212)/sum(g212_A%g212_A));
        double res222 = sqrt(sum(res_g222%res_g222)/sum(g222_A%g222_A));
        std::cout << "nx = " << nx << "\tny = " << ny << '\n';
        std::cout << "res_e_11 = " << res_11 << "\tres_e_12 = " << res_12 << "\tres_e_22 = " << res_22 << '\n';
        std::cout << "res_e11  = " << res11  << "\tres_e12  = " << res12  << "\tres_e22  = " << res22  << '\n';
        std::cout << "res_g111 = " << res111 << "\tres_g112 = " << res112 << "\tres_g122 = " << res122 << '\n';
        std::cout << "res_g211 = " << res211 << "\tres_g212 = " << res212 << "\tres_g222 = " << res222 << "\n\n";
        
        file << nx << ' ' << ny << ' '
                << res_11 << ' ' << res_12 << ' ' << res_22 << ' '
                << res11  << ' ' << res12  << ' ' << res22  << ' '
                << res111 << ' ' << res112 << ' ' << res122 << ' '
                << res211 << ' ' << res212 << ' ' << res222 << '\n';

        std::ofstream file_xys ("plot/Data/xys");
        for (size_t i = 0; i < nx; i++, file_xys<<'\n')
            for (size_t j = 0; j < ny; j++, file_xys<<'\n')
                file_xys  << xs(i, j)  << ' ' << ys(i, j)  << " 0";
        file_xys.close();
    }
    file.close();
}