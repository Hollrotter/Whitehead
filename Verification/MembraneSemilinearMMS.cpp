#include "Membrane.hpp"
#include "interpolation.hpp"
#include <ginac/ginac.h>

using namespace GiNaC;

int main()
{
    symbol x1("x1"), x2("x2");
    
    size_t nx = 20;
    size_t ny = 20;
    size_t samples = 1;

    const double pi = arma::datum::pi;

    const double r1_N = 1;
    const double r1_S = 1;
    const double r1_W = 1;
    const double r1_E = 1;
    const double r2_N = 0;
    const double r2_S = 0;
    const double r2_W = 0;
    const double r2_E = 0;

    // Geometry
    
    // ex x = x1 + 0.1*sin(pi*x2)*cos(pi/8*x1);
    // ex y = x2 + 0.1*sin(pi*x1)*cos(pi/8*x2);
    double a1 = 0.125;
    double a2 = 0.125;

    ex x = x1*(1 + a1*sin(pi*x2));
    ex y = x2*(1 + a2*sin(pi*x1));

    // Displacement
    // ex z   = (1 - pow(x1, 2))*(1 - pow(x2, 2));
    ex z   = 1 - pow(x, 2) - pow(y, 2);
    // ex z   = (1 - pow(x1, 2) - pow(x2, 2))/4;
    
    ex kxx = 1.0 + sin(pi*x) * cos(pi*y);
    ex kxy = 0.1 * cos(pi*x) * cos(pi*y);
    ex kyy = 1.0 + cos(pi*x) * sin(pi*y);
    
    ex dxdx1 = diff(x, x1);
    ex dxdx2 = diff(x, x2);
    ex dydx1 = diff(y, x1);
    ex dydx2 = diff(y, x2);
    
    ex det2 = pow(dxdx1*dydx2 - dxdx2*dydx1, 2);

    ex k11 = (dydx2*dydx2*kxx -              2*dydx2*dxdx2 *kxy + dxdx2*dxdx2*kyy)/det2;
    ex k12 =-(dydx2*dydx1*kxx - (dxdx1*dydx2 + dydx1*dxdx2)*kxy + dxdx2*dxdx1*kyy)/det2;
    ex k22 = (dydx1*dydx1*kxx -              2*dydx1*dxdx1 *kxy + dxdx1*dxdx1*kyy)/det2;

    ex e_11 = dxdx1*dxdx1 + dydx1*dydx1;
    ex e_12 = dxdx1*dxdx2 + dydx1*dydx2;
    ex e_22 = dxdx2*dxdx2 + dydx2*dydx2;
    
    ex e = e_11*e_22 - e_12*e_12;
    
    ex e11 = e_22/e;
    ex e12 =-e_12/e;
    ex e22 = e_11/e;
    
    ex de_11dx1 = diff(e_11, x1);
    ex de_12dx1 = diff(e_12, x1);
    ex de_22dx1 = diff(e_22, x1);
    ex de_11dx2 = diff(e_11, x2);
    ex de_12dx2 = diff(e_12, x2);
    ex de_22dx2 = diff(e_22, x2);
    
    ex g111 = e11/2*de_11dx1 + e12*(de_12dx1 - de_11dx2/2);
    ex g112 = e11/2*de_11dx2 + e12/2*de_22dx1;
    ex g122 = e11*(de_12dx2 - de_22dx1/2) + e12/2*de_22dx2;
    ex g211 = e12/2*de_11dx1 + e22*(de_12dx1 - de_11dx2/2);
    ex g212 = e12/2*de_11dx2 + e22/2*de_22dx1;
    ex g222 = e12*(de_12dx2 - de_22dx1/2) + e22/2*de_22dx2;
    
    // PDE
    ex z_1  = diff(z, x1);
    ex z_2  = diff(z, x2);
    ex z_11 = diff(z_1, x1);
    ex z_12 = diff(z_1, x2);
    ex z_22 = diff(z_2, x2);
    ex z__11 = z_11 - g111*z_1 - g211*z_2;
    ex z__12 = z_12 - g112*z_1 - g212*z_2;
    ex z__22 = z_22 - g122*z_1 - g222*z_2;
    
    ex sae = sqrt(1 + e11*z_1*z_1 + 2*e12*z_1*z_2 + e22*z_2*z_2);
    ex p =-(z__11*k11 + 2*z__12*k12 + z__22*k22)/sae;

    ex det = dxdx1*dydx2 - dxdx2*dydx1;
    ex z_x =( dydx2*z_1 - dydx1*z_2)/det;
    ex z_y =(-dxdx2*z_1 + dxdx1*z_2)/det;

    Membrane m;
    
    std::ofstream file("plot/Data/SemilinearMMS");

    for (size_t k = 0; k < samples; k++, nx+=2, ny+=2)
    {
        arma::vec x1s = Chebyshev::gaussLobatto(nx);
        arma::vec x2s = Chebyshev::gaussLobatto(ny);

        arma::mat ys(nx, ny), xs(nx, ny), ps(nx, ny), kxxs(nx, ny), kxys(nx, ny), kyys(nx, ny);
        arma::vec z_A(nx*ny), zx_A(nx*ny), zy_A(nx*ny);

        exmap mp;
        for (size_t i = 0; i < nx; i++)
        {
            std::cout << "i = " << i << "\n\tj = ";
            mp[x1] = x1s(i);
            for (size_t j = 0; j < ny; j++)
            {
                std::cout << j << ' ';
                mp[x2] = x2s(j);
                xs(i, j)     = ex_to<numeric>(evalf(x  .subs(mp))).to_double();
                ys(i, j)     = ex_to<numeric>(evalf(y  .subs(mp))).to_double();
                ps(i, j)     = ex_to<numeric>(evalf(p  .subs(mp))).to_double();
                kxxs(i, j)   = ex_to<numeric>(evalf(kxx.subs(mp))).to_double();
                kxys(i, j)   = ex_to<numeric>(evalf(kxy.subs(mp))).to_double();
                kyys(i, j)   = ex_to<numeric>(evalf(kyy.subs(mp))).to_double();
                z_A(i+j*nx)  = ex_to<numeric>(evalf(z  .subs(mp))).to_double();
                zx_A(i+j*nx) = ex_to<numeric>(evalf(z_x.subs(mp))).to_double();
                zy_A(i+j*nx) = ex_to<numeric>(evalf(z_y.subs(mp))).to_double();
            }
            std::cout << '\n';
        }
        Lagrange::CurveInterpolant chi1(xs.col(0),        ys.col(0));
        Lagrange::CurveInterpolant chi2(xs.row(0).t(),    ys.row(0).t());
        Lagrange::CurveInterpolant chi3(xs.col(ny-1),     ys.col(ny-1));
        Lagrange::CurveInterpolant chi4(xs.row(nx-1).t(), ys.row(nx-1).t());
        
        auto [xss, yss] = Lagrange::TransfiniteQuadMap({&chi1, &chi2, &chi3, &chi4});

        arma::mat T = ChebyshevInterpolation(x1s, x2s, xs, ys, xss, yss);
        z_A  = T*z_A;
        zx_A = T*zx_A;
        zy_A = T*zy_A;
        arma::mat kxxss = arma::reshape(T*vectorise(kxxs), nx, ny);
        arma::mat kxyss = arma::reshape(T*vectorise(kxys), nx, ny);
        arma::mat kyyss = arma::reshape(T*vectorise(kyys), nx, ny);
        arma::mat pss   = arma::reshape(T*vectorise(ps),   nx, ny);

        m({&chi1, &chi2, &chi3, &chi4});
        arma::mat e11ss = m.e11();
        arma::mat e12ss = m.e12();
        arma::mat e22ss = m.e22();
        arma::mat dxdx1ss = m.J11();
        arma::mat dxdx2ss = m.J12();
        arma::mat dydx1ss = m.J21();
        arma::mat dydx2ss = m.J22();
        arma::mat z_1_A = dxdx1ss%reshape(zx_A, nx, ny) + dydx1ss%reshape(zy_A, nx, ny);
        arma::mat z_2_A = dxdx2ss%reshape(zx_A, nx, ny) + dydx2ss%reshape(zy_A, nx, ny);

        // Boundary Conditions
        // West and East
        
        arma::mat h_1s1 = sqrt(e11ss);
        arma::mat h_1s2 = e12ss/sqrt(e11ss);
        
        arma::mat z_1s_WE = h_1s1%z_1_A + h_1s2%z_2_A;

        // North and South
        
        arma::mat h_2s1 = e12ss/sqrt(e22ss);
        arma::mat h_2s2 = sqrt(e22ss);

        arma::mat z_2s_NS = h_2s1%z_1_A + h_2s2%z_2_A;

        arma::vec z_N(nx), z_S(nx), z_W(ny), z_E(ny);

        for (size_t i = 0; i < nx; i++)
        {
            z_S(i)  = r1_S*z_A(i)           + r2_S*z_2s_NS(i, 0);
            z_N(i)  = r1_N*z_A(i+nx*(ny-1)) + r2_N*z_2s_NS(i, ny-1);
        }
        for (size_t j = 0; j < ny; j++)
        {
            z_W(j)  = r1_W*z_A(nx*j)      + r2_W*z_1s_WE(0,    j);
            z_E(j)  = r1_E*z_A(nx-1+nx*j) + r2_E*z_1s_WE(nx-1, j);
        }

        m.boundary(Field::z, Direction::N, BC::Robin, r1_N, r2_N, z_N);
        m.boundary(Field::z, Direction::S, BC::Robin, r1_S, r2_S, z_S);
        m.boundary(Field::z, Direction::W, BC::Robin, r1_W, r2_W, z_W);
        m.boundary(Field::z, Direction::E, BC::Robin, r1_E, r2_E, z_E);

        m.load(pss);
        m.inPlaneKartesian(kxxss, kxyss, kyyss);

        m.semilinear();

        arma::vec z_num = vectorise(arma::mat(m(Field::z)));
        
        m.output("plot/Data/z", Field::z);
        std::ofstream file_z("plot/Data/z_A");
        std::ofstream file_zres("plot/Data/z_res");
        for (size_t i = 0; i < nx; i++, file_z<<'\n', file_zres<<'\n')
            for (size_t j = 0; j < ny; j++, file_z<<'\n', file_zres<<'\n')
            {
                file_z    << xss(i, j) << ' ' << yss(i, j) << ' ' << z_A(i+nx*j);
                file_zres << xss(i, j) << ' ' << yss(i, j) << ' ' << z_num(i+nx*j)-z_A(i+nx*j);
            }
        file_z.close();
        file_zres.close();

        arma::vec res_z = z_num - z_A;
        double res = sqrt(sum(res_z%res_z)/sum(z_A%z_A));
        std::cout << "nx = " << nx << " ny = " << ny << " res_z = " << res << '\n';

        file << nx << ' ' << ny << ' ' << res << '\n';
    }
    file.close();
}