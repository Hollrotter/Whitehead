#include "Membrane.hpp"
#include "interpolation.hpp"
#include <ginac/ginac.h>

using namespace GiNaC;

int main()
{
    symbol x1("x1"), x2("x2");
    
    size_t nx = 21;
    size_t ny = 21;
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

    double a1 = 0.15;
    double a2 = 0.1;

    ex x = x1*(1 + a1*sin(pi*x2));
    ex y = x2*(1 + a2*sin(pi*x1));

    // ex x = x1 + 0.15*sin(pi*x2);
    // ex y = x2 + 0.15*sin(pi*x1);

    // ex x = (x1+2)*cos(x2);
    // ex y = (x1+2)*sin(x2);

    // Displacement
    // ex z = 1 - pow(x, 2) - pow(y, 2);
    ex z = cos(x)*cos(y);

    ex nxx = 1.0 + sin(pi*x) * cos(pi*y);
    ex nxy = 0.1 * cos(pi*x) * cos(pi*y);
    ex nyy = 1.0 + cos(pi*x) * sin(pi*y);

    ex dxdx1 = diff(x, x1);
    ex dxdx2 = diff(x, x2);
    ex dydx1 = diff(y, x1);
    ex dydx2 = diff(y, x2);

    ex det2 = pow(dxdx1*dydx2 - dxdx2*dydx1, 2);

    ex n11 = (dydx2*dydx2*nxx -              2*dxdx2*dydx2 *nxy + dxdx2*dxdx2*nyy)/det2;
    ex n12 =-(dydx1*dydx2*nxx - (dydx1*dxdx2 + dxdx1*dydx2)*nxy + dxdx1*dxdx2*nyy)/det2;
    ex n22 = (dydx1*dydx1*nxx -              2*dxdx1*dydx1 *nxy + dxdx1*dxdx1*nyy)/det2;

    ex e_11 = dxdx1*dxdx1 + dydx1*dydx1;
    ex e_12 = dxdx1*dxdx2 + dydx1*dydx2;
    ex e_22 = dxdx2*dxdx2 + dydx2*dydx2;
    
    ex e = e_11*e_22 - pow(e_12, 2);

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
    
    // std::cout << g111 << '\n';
    // std::cout << g112 << '\n';
    // std::cout << g122 << '\n';
    // std::cout << g211 << '\n';
    // std::cout << g212 << '\n';
    // std::cout << g222 << '\n';

    // PDE
    ex z_1  = diff(z, x1);
    ex z_2  = diff(z, x2);
    ex z_11 = diff(z_1, x1);
    ex z_12 = diff(z_1, x2);
    ex z_22 = diff(z_2, x2);
    ex z__11 = z_11 - g111*z_1 - g211*z_2;
    ex z__12 = z_12 - g112*z_1 - g212*z_2;
    ex z__22 = z_22 - g122*z_1 - g222*z_2;
    
    ex p =-(z__11*n11 + 2*z__12*n12 + z__22*n22);
    
    ex det = dxdx1*dydx2 - dxdx2*dydx1;
    ex z_x =( dydx2*z_1 - dydx1*z_2)/det;
    ex z_y =(-dxdx2*z_1 + dxdx1*z_2)/det;

    Membrane m;

    std::ofstream file("plot/Data/MembraneLinearMMS");

    for (size_t k = 0; k < samples; k++, nx+=2, ny+=2)
    {
        arma::vec x1s = Chebyshev::gaussLobatto(nx);
        arma::vec x2s = Chebyshev::gaussLobatto(ny);

        arma::mat ys(nx, ny), xs(nx, ny), ps(nx, ny), nxxs(nx, ny), nxys(nx, ny), nyys(nx, ny);
        arma::vec z_A(nx*ny), zx_A(nx*ny), zy_A(nx*ny);

        exmap mp;
        for (size_t i = 0; i < nx; i++)
        {
            std::cout << "i = " << i << '\n';
            mp[x1] = x1s(i);
            for (size_t j = 0; j < ny; j++)
            {
                mp[x2] = x2s(j);
                xs(i, j)     = ex_to<numeric>(evalf(x  .subs(mp))).to_double();
                ys(i, j)     = ex_to<numeric>(evalf(y  .subs(mp))).to_double();
                ps(i, j)     = ex_to<numeric>(evalf(p  .subs(mp))).to_double();
                nxxs(i, j)   = ex_to<numeric>(evalf(nxx.subs(mp))).to_double();
                nxys(i, j)   = ex_to<numeric>(evalf(nxy.subs(mp))).to_double();
                nyys(i, j)   = ex_to<numeric>(evalf(nyy.subs(mp))).to_double();
                z_A(i+j*nx)  = ex_to<numeric>(evalf(z  .subs(mp))).to_double();
                zx_A(i+j*nx) = ex_to<numeric>(evalf(z_x.subs(mp))).to_double();
                zy_A(i+j*nx) = ex_to<numeric>(evalf(z_y.subs(mp))).to_double();
            }
        }

        Lagrange::CurveInterpolant chi1(xs.col(0),        ys.col(0));
        Lagrange::CurveInterpolant chi2(xs.row(nx-1).t(), ys.row(nx-1).t());
        Lagrange::CurveInterpolant chi3(xs.col(ny-1),     ys.col(ny-1));
        Lagrange::CurveInterpolant chi4(xs.row(0).t(),    ys.row(0).t());
        
        auto [xss, yss] = Lagrange::TransfiniteQuadMap({&chi1, &chi2, &chi3, &chi4});

        arma::mat T = ChebyshevInterpolation(x1s, x2s, xs, ys, xss, yss);
        z_A  = T*z_A;
        zx_A = T*zx_A;
        zy_A = T*zy_A;
        arma::mat nxxss = arma::reshape(T*vectorise(nxxs), nx, ny);
        arma::mat nxyss = arma::reshape(T*vectorise(nxys), nx, ny);
        arma::mat nyyss = arma::reshape(T*vectorise(nyys), nx, ny);
        arma::mat pss   = arma::reshape(T*vectorise(ps),   nx, ny);

        m({&chi1, &chi2, &chi3, &chi4});

        // m.g111().print();
        // m.g112().print();
        // m.g122().print();
        // m.g211().print();
        // m.g212().print();
        // m.g222().print();
        arma::mat z_1_A = m.J11()%reshape(zx_A, nx, ny) + m.J21()%reshape(zy_A, nx, ny);
        arma::mat z_2_A = m.J12()%reshape(zx_A, nx, ny) + m.J22()%reshape(zy_A, nx, ny);

        // Boundary Conditions
        arma::mat e11ss = m.e11();
        arma::mat e12ss = m.e12();
        arma::mat e22ss = m.e22();

        // West and East
        arma::mat h_1s1 = sqrt(e11ss);
        arma::mat h_1s2 = e12ss/sqrt(e11ss);
        
        arma::mat z_1s_WE = h_1s1%z_1_A + h_1s2%z_2_A;

        // North and South
        arma::mat h_2s1 = e12ss/sqrt(e22ss);
        arma::mat h_2s2 = sqrt(e22ss);

        arma::mat z_2s_NS = h_2s1%z_1_A + h_2s2%z_2_A;

        arma::vec z_W(ny), z_E(ny), z_S(nx), z_N(nx);

        for (size_t i = 0; i < nx; i++)
        {
            z_S(i) = r1_S*z_A(i)           + r2_S*z_2s_NS(i,    0);
            z_N(i) = r1_N*z_A(i+nx*(ny-1)) + r2_N*z_2s_NS(i, ny-1);
        }
        for (size_t j = 0; j < ny; j++)
        {
            z_W(j) = r1_W*z_A(nx*j)       + r2_W*z_1s_WE(   0, j);
            z_E(j) = r1_E*z_A(nx*(1+j)-1) + r2_E*z_1s_WE(nx-1, j);
        }

        m.boundary(Field::z, Direction::N, BC::Robin, r1_N, r2_N, z_N);
        m.boundary(Field::z, Direction::S, BC::Robin, r1_S, r2_S, z_S);
        m.boundary(Field::z, Direction::W, BC::Robin, r1_W, r2_W, z_W);
        m.boundary(Field::z, Direction::E, BC::Robin, r1_E, r2_E, z_E);

        m.load(pss);
        m.inPlaneKartesian(nxxss, nxyss, nyyss);

        m.linear();

        arma::vec z_num = vectorise(arma::mat(m(Field::z)));
        
        // m.output("plot/Data/z", Field::z);
        // std::ofstream file_z("plot/Data/z_A");
        // std::ofstream file_zres("plot/Data/z_res");
        // for (size_t i = 0; i < nx; i++, file_z<<'\n', file_zres<<'\n')
        //     for (size_t j = 0; j < ny; j++, file_z<<'\n', file_zres<<'\n')
        //     {
        //         file_z    << xss(i, j) << ' ' << yss(i, j) << ' ' << z_A(i+nx*j);
        //         file_zres << xss(i, j) << ' ' << yss(i, j) << ' ' << z_num(i+nx*j)-z_A(i+nx*j);
        //     }
        // file_z.close();
        // file_zres.close();

        arma::vec res_z = z_num - z_A;
        double res = sqrt(sum(res_z%res_z)/sum(z_A%z_A));
        std::cout << "nx = " << nx << " ny = " << ny << " res_z = " << res << '\n';

        file << nx << ' ' << ny << ' ' << res << '\n';

        std::ofstream file_xys ("plot/Data/xys");
        std::ofstream file_xyss("plot/Data/xyss");

        for (size_t i = 0; i < nx; i++, file_xys<<'\n', file_xyss<<'\n')
            for (size_t j = 0; j < ny; j++, file_xys<<'\n', file_xyss<<'\n')
            {
                file_xys  << xs(i, j)  << ' ' << ys(i, j)  << " 0";
                file_xyss << xss(i, j) << ' ' << yss(i, j) << " 0";
            }
        file_xys.close();
        file_xyss.close();
    }

    file.close();
}