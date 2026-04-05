#include "Membrane.hpp"
#include "interpolation.hpp"
#include <ginac/ginac.h>

using namespace GiNaC;

int main()
{
    symbol x1("x1"), x2("x2");
    
    size_t nx = 10;
    size_t ny = 10;
    size_t samples = 1;

    const double nu = 0.3;
    const double Et = 1000;
    const double D = Et/(1 - nu*nu);
    const double pi = arma::datum::pi;

    // const double r1_N = 1;
    // const double r1_S = 1;
    // const double r1_W = 1;
    // const double r1_E = 1;
    // const double r2_N = 1;
    // const double r2_S = 1;
    // const double r2_W = 1;
    // const double r2_E = 1;

    // Geometry
    // ex x = x1 + 0.1*sin(pi*x2)*cos(pi/8*x1);
    // ex y = x2 + 0.1*sin(pi*x1)*cos(pi/8*x2);
    double a1 = 0.125;
    double a2 = 0.125;

    ex x = x1*(1 + a1*sin(pi*x2));
    ex y = x2*(1 + a2*sin(pi*x1));

    // Displacements
    ex v_x = x1;//sin(pi*x);//*cos(pi/8*y);// + 0.1*sin(pi*x2)*cos(pi/8*x2);
    ex v_y = x2;//sin(pi*y);//*cos(pi/8*x);// + 0.1*sin(pi*x1)*cos(pi/8*x1);

    // ex v_x = sin(pi*x1)*cos(pi/8*x2);
    // ex v_y = sin(pi*x2)*cos(pi/8*x1);
    
    ex dxdx1 = diff(x, x1);
    ex dxdx2 = diff(x, x2);
    ex dydx1 = diff(y, x1);
    ex dydx2 = diff(y, x2);

    ex v_1 = dxdx1*v_x + dydx1*v_y;
    ex v_2 = dxdx2*v_x + dydx2*v_y;

    ex v_1_1  = diff(v_1,   x1);
    ex v_1_2  = diff(v_1,   x2);
    ex v_2_1  = diff(v_2,   x1);
    ex v_2_2  = diff(v_2,   x2);
    
    ex v_1_11 = diff(v_1_1, x1);
    ex v_1_12 = diff(v_1_1, x2);
    ex v_1_22 = diff(v_1_2, x2);
    ex v_2_11 = diff(v_2_1, x1);
    ex v_2_12 = diff(v_2_1, x2);
    ex v_2_22 = diff(v_2_2, x2);
    
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
    
    ex g111_1 = diff(g111, x1);
    ex g112_1 = diff(g112, x1);
    ex g122_1 = diff(g122, x1);
    ex g211_1 = diff(g211, x1);
    ex g212_1 = diff(g212, x1);
    ex g222_1 = diff(g222, x1);
    
    ex g111_2 = diff(g111, x2);
    ex g112_2 = diff(g112, x2);
    ex g122_2 = diff(g122, x2);
    ex g211_2 = diff(g211, x2);
    ex g212_2 = diff(g212, x2);
    ex g222_2 = diff(g222, x2);
    
    // PDEs
    ex v_1__1 = v_1_1 - g111*v_1 - g211*v_2;
    ex v_1__2 = v_1_2 - g112*v_1 - g212*v_2;
    ex v_2__1 = v_2_1 - g112*v_1 - g212*v_2;
    ex v_2__2 = v_2_2 - g122*v_1 - g222*v_2;
    
    ex v_1__11 = v_1_11 - 3*g111*v_1_1 - g211*v_1_2 - 2*g211*v_2_1
               + (2*g111*g111 + 2*g211*g112 - g111_1)*v_1
               + (2*g111*g211 + 2*g211*g212 - g211_1)*v_2;
    ex v_1__12 = v_1_12 - 2*g112*v_1_1 - g212*v_2_1 - g211*v_2_2 - (g111 + g212)*v_1_2
               + (2*g112*g111 + 2*g212*g112 - g111_2)*v_1
               + (2*g112*g211 + 2*g212*g212 - g211_2)*v_2;
    ex v_1__22 = v_1_22 - g122*v_1_1 - (2*g112 + g222)*v_1_2 - 2*g212*v_2_2
               + (g122*g111 + g222*g112 + g112*g112 + g212*g122 - g112_2)*v_1
               + (g122*g211 + g112*g212 + 2*g212*g222 - g212_2)*v_2;
    ex v_2__11 = v_2_11 - 2*g112*v_1_1 - (2*g212 + g111)*v_2_1 - g211*v_2_2
               + (2*g111*g112 + g211*g122 + g212*g112 - g112_1)*v_1
               + (g111*g212 + g211*g222 + g112*g211 + g212*g212 - g212_1)*v_2;
    ex v_2__12 = v_2_12 - g122*v_1_1 - g112*v_1_2 - 2*g212*v_2_2 - (g112 + g222)*v_2_1
               + (2*g112*g112 + 2*g212*g122 - g122_1)*v_1
               + (2*g112*g212 + 2*g212*g222 - g222_1)*v_2;
    ex v_2__22 = v_2_22 - 2*g122*v_1_2 - g122*v_2_1 - 3*g222*v_2_2
               + (2*g122*g112 + 2*g222*g122 - g122_2)*v_1
               + (2*g122*g212 + 2*g222*g222 - g222_2)*v_2;
    
    ex g_11 = v_1__1;
    ex g_12 = (v_1__2 + v_2__1)/2;
    ex g_22 = v_2__2;
    
    ex g_11_1 = v_1__11;
    ex g_11_2 = v_1__12;
    ex g_12_1 = (v_1__12 + v_2__11)/2;
    ex g_12_2 = (v_1__22 + v_2__12)/2;
    ex g_22_1 = v_2__12;
    ex g_22_2 = v_2__22;
    
    ex H1111 = e11*e11;
    ex H1112 = e11*e12;
    ex H1122 = e12*e12 + nu/e;
    ex H1221 = (e11*e22 + e12*e12 - nu/e)/2;
    ex H1222 = e12*e22;
    ex H2222 = e22*e22;
    
    ex n11 = D*(H1111*g_11 + 2*H1112*g_12 + H1122*g_22);
    ex n12 = D*(H1112*g_11 + 2*H1221*g_12 + H1222*g_22);
    ex n22 = D*(H1122*g_11 + 2*H1222*g_12 + H2222*g_22);
    
    ex n11__1 = D*(H1111*g_11_1 + 2*H1112*g_12_1 + H1122*g_22_1);
    ex n12__1 = D*(H1112*g_11_1 + 2*H1221*g_12_1 + H1222*g_22_1);
    ex n12__2 = D*(H1112*g_11_2 + 2*H1221*g_12_2 + H1222*g_22_2);
    ex n22__2 = D*(H1122*g_11_2 + 2*H1222*g_12_2 + H2222*g_22_2);
    
    ex p1 =-n11__1 - n12__2;
    ex p2 =-n12__1 - n22__2;

    ex px = dxdx1*p1 + dxdx2*p2;
    ex py = dydx1*p1 + dydx2*p2;

    ex nxx = dxdx1*dxdx1*n11 +               2*dxdx1*dxdx2*n12 + dxdx2*dxdx2*n22;
    ex nxy = dxdx1*dydx1*n11 + (dxdx1*dydx2 + dxdx2*dydx1)*n12 + dxdx2*dydx2*n22;
    ex nyy = dydx1*dydx1*n11 +               2*dydx1*dydx2*n12 + dydx2*dydx2*n22;

    std::ofstream file("plot/Data/MembranePlaneStrainMMS");

    for (size_t k = 0; k < samples; k++, nx+=2, ny+=2)
    {
        arma::vec x1s = Chebyshev::gaussLobatto(nx);
        arma::vec x2s = Chebyshev::gaussLobatto(ny);

        arma::mat xs(nx, ny), ys(nx, ny);
        arma::vec vx_A(nx*ny), vy_A(nx*ny), pxs(nx*ny), pys(nx*ny);
        arma::vec nxx_A(nx*ny), nxy_A(nx*ny), nyy_A(nx*ny);

        exmap mp;
        for (size_t i = 0; i < nx; i++)
        {
            std::cout << "i = " << i << '\n';
            mp[x1] = x1s(i);
            std::cout << '\t' << "j = ";
            for (size_t j = 0; j < ny; j++)
            {
                std::cout << j << ' ';
                mp[x2] = x2s(j);
                xs(i, j)      = ex_to<numeric>(evalf(x  .subs(mp))).to_double();
                ys(i, j)      = ex_to<numeric>(evalf(y  .subs(mp))).to_double();
                pxs(i+j*nx)   = ex_to<numeric>(evalf(px .subs(mp))).to_double();
                pys(i+j*nx)   = ex_to<numeric>(evalf(py .subs(mp))).to_double();
                vx_A(i+j*nx)  = ex_to<numeric>(evalf(v_x.subs(mp))).to_double();
                vy_A(i+j*nx)  = ex_to<numeric>(evalf(v_y.subs(mp))).to_double();
                nxx_A(i+j*nx) = ex_to<numeric>(evalf(nxx.subs(mp))).to_double();
                nxy_A(i+j*nx) = ex_to<numeric>(evalf(nxy.subs(mp))).to_double();
                nyy_A(i+j*nx) = ex_to<numeric>(evalf(nyy.subs(mp))).to_double();
            }
            std::cout << '\n';
        }

        Lagrange::CurveInterpolant chi1(xs.col(0),        ys.col(0));
        Lagrange::CurveInterpolant chi2(xs.row(nx-1).t(), ys.row(nx-1).t());
        Lagrange::CurveInterpolant chi3(xs.col(ny-1),     ys.col(ny-1));
        Lagrange::CurveInterpolant chi4(xs.row(0).t(),    ys.row(0).t());
        
        auto [xss, yss] = Lagrange::TransfiniteQuadMap({&chi1, &chi2, &chi3, &chi4});
        arma::mat T = ChebyshevInterpolation(x1s, x2s, xs, ys, xss, yss);
        vx_A  = T * vx_A;
        vy_A  = T * vy_A;
        nxx_A = T * nxx_A;
        nxy_A = T * nxy_A;
        nyy_A = T * nyy_A;

        arma::vec pxss = T * pxs;
        arma::vec pyss = T * pys;

        // Membrane m(xss, yss);
        Membrane m({&chi1, &chi2, &chi3, &chi4});
        m.youngsModulus(Et);
        m.poissonsRatio(nu);
        
        arma::mat e_11ss = m.e_11();
        arma::mat e_12ss = m.e_12();
        arma::mat e_22ss = m.e_22();
        arma::mat e11ss  = m.e11();
        arma::mat e12ss  = m.e12();
        arma::mat e22ss  = m.e22();

        arma::mat dxdx1ss = m.J11();
        arma::mat dxdx2ss = m.J12();
        arma::mat dydx1ss = m.J21();
        arma::mat dydx2ss = m.J22();

        arma::mat det2_ss = (dxdx1ss%dydx2ss - dydx1ss%dxdx2ss)%(dxdx1ss%dydx2ss - dydx1ss%dxdx2ss);
        
        arma::mat v1_A = dxdx1ss%reshape(vx_A, nx, ny) + dydx1ss%reshape(vy_A, nx, ny);
        arma::mat v2_A = dxdx2ss%reshape(vx_A, nx, ny) + dydx2ss%reshape(vy_A, nx, ny);

        arma::mat n11_A = (dydx2ss%dydx2ss%reshape(nxx_A, nx, ny) -                   2*dydx2ss%dxdx2ss%reshape(nxy_A, nx, ny) + dxdx2ss%dxdx2ss%reshape(nyy_A, nx, ny))/det2_ss;
        arma::mat n12_A =-(dydx2ss%dydx1ss%reshape(nxx_A, nx, ny) - (dydx2ss%dxdx1ss + dydx1ss%dxdx2ss)%reshape(nxy_A, nx, ny) + dxdx2ss%dxdx1ss%reshape(nyy_A, nx, ny))/det2_ss;
        arma::mat n22_A = (dydx1ss%dydx1ss%reshape(nxx_A, nx, ny) -                   2*dydx1ss%dxdx1ss%reshape(nxy_A, nx, ny) + dxdx1ss%dxdx1ss%reshape(nyy_A, nx, ny))/det2_ss;
        
        // Boundary Conditions
            
        arma::vec v1_N(nx), v1_S(nx), v1_W(ny), v1_E(ny);
        arma::vec v2_N(nx), v2_S(nx), v2_W(ny), v2_E(ny);
        arma::vec n11_W(ny), n11_E(ny);
        arma::vec n12_W(ny), n12_E(ny);
        arma::vec n21_N(nx), n21_S(nx);
        arma::vec n22_N(nx), n22_S(nx);

        // West and East
        arma::mat h_11s = sqrt(1/e11ss);
        arma::mat h_12s = e_12ss/sqrt(e_22ss);
        arma::mat h_22s = sqrt(e_22ss);
        
        arma::mat h_1s1 = 1/h_11s;
        arma::mat h_1s2 = h_11s%e12ss;
        arma::mat h_2s2 = 1/h_22s;
        
        arma::mat v_1s_WE = h_1s1%v1_A + h_1s2%v2_A;
        arma::mat v_2s_WE = h_2s2%v2_A;
        
        // arma::mat v_1s__1s = h_1s1*h_1s1*v_1__1 + h_1s1*h_1s2*(v_1__2 + v_2__1) + h_1s2*h_1s2*v_2__2;
        // arma::mat v_2s__1s = h_1s1*h_2s2*v_2__1 + h_2s2*h_1s2*v_2__2;

        arma::mat n_1s1s_WE = h_11s%h_11s%n11_A;
        arma::mat n_1s2s_WE = h_11s%h_12s%n11_A + h_11s%h_22s%n12_A;
        
        for (size_t j = 0; j < ny; j++)
        {
            v1_W(j)  = v_1s_WE(0, j);
            v2_W(j)  = v_2s_WE(0, j);
            v1_E(j)  = v_1s_WE(nx-1, j);
            v2_E(j)  = v_2s_WE(nx-1, j);
            n11_W(j) = n_1s1s_WE(0, j);
            n12_W(j) = n_1s2s_WE(0, j);
            n11_E(j) = n_1s1s_WE(nx-1, j);
            n12_E(j) = n_1s2s_WE(nx-1, j);
        }

        // North and South
        h_11s = sqrt(e_11ss);
        arma::mat h_21s = e_12ss/sqrt(e_11ss);
        h_22s = sqrt(1/e22ss);
        
        h_1s1 = 1/h_11s;
        arma::mat h_2s1 = h_22s%e12ss;
        h_2s2 = 1/h_22s;
        
        arma::mat v_1s_NS = h_1s1%v1_A;
        arma::mat v_2s_NS = h_2s1%v1_A + h_2s2%v2_A;
        
        // arma::mat v_1s__2s = h_1s1*h_2s1*v_1__1 + h_1s1*h_2s2*v_1__2;
        // arma::mat v_2s__2s = h_2s1*h_2s1*v_1__1 + h_2s1*h_2s2*(v_1__2 + v_2__1) + h_2s2*h_2s2*v_2__2;
        
        arma::mat n_2s2s_NS = h_22s%h_22s%n22_A;
        arma::mat n_2s1s_NS = h_11s%h_22s%n12_A + h_21s%h_22s%n22_A;

        for (size_t i = 0; i < nx; i++)
        {
            v1_S(i)  = v_1s_NS(i, 0);
            v2_S(i)  = v_2s_NS(i, 0);
            v1_N(i)  = v_1s_NS(i, ny-1);
            v2_N(i)  = v_2s_NS(i, ny-1);
            n22_S(i) = n_2s2s_NS(i, 0);
            n21_S(i) = n_2s1s_NS(i, 0);
            n22_N(i) = n_2s2s_NS(i, ny-1);
            n21_N(i) = n_2s1s_NS(i, ny-1);
        }
        // m.boundary(Field::v1, Direction::N, BC::Dirichlet, v1_N);
        m.boundary(Field::v1, Direction::S, BC::Dirichlet, v1_S);
        m.boundary(Field::v1, Direction::W, BC::Dirichlet, v1_W);
        m.boundary(Field::v1, Direction::E, BC::Dirichlet, v1_E);

        m.boundary(Field::n12, Direction::N, BC::Dirichlet, n21_N);
        // m.boundary(Field::n12, Direction::S, BC::Dirichlet, n21_S);
        // m.boundary(Field::n11, Direction::W, BC::Dirichlet, n11_W);
        // m.boundary(Field::n11, Direction::E, BC::Dirichlet, n11_E);

        // m.boundary(Field::v2, Direction::N, BC::Dirichlet, v2_N);
        m.boundary(Field::v2, Direction::S, BC::Dirichlet, v2_S);
        m.boundary(Field::v2, Direction::W, BC::Dirichlet, v2_W);
        m.boundary(Field::v2, Direction::E, BC::Dirichlet, v2_E);

        m.boundary(Field::n22, Direction::N, BC::Dirichlet, n22_N);
        // m.boundary(Field::n22, Direction::S, BC::Dirichlet, n22_S);
        // m.boundary(Field::n12, Direction::W, BC::Dirichlet, n12_W);
        // m.boundary(Field::n12, Direction::E, BC::Dirichlet, n12_E);

        m.inPlaneKartesian(pxss, pyss);

        m.planeStrain();
        
        auto [vxss, vyss] = m.kartesianDisplacements();
        arma::vec res_vx = vectorise(vxss) - vx_A;
        arma::vec res_vy = vectorise(vyss) - vy_A;

        std::ofstream file_vx("plot/Data/vx");
        std::ofstream file_vy("plot/Data/vy");
        std::ofstream file_vx_A("plot/Data/vx_A");
        std::ofstream file_vy_A("plot/Data/vy_A");
        std::ofstream file_vxres("plot/Data/vx_res");
        std::ofstream file_vyres("plot/Data/vy_res");
        for (size_t i = 0; i < nx; i++, file_vx<<'\n', file_vx_A<<'\n' , file_vxres<<'\n', file_vy<<'\n', file_vy_A<<'\n' , file_vyres<<'\n')
            for (size_t j = 0; j < ny; j++, file_vx<<'\n', file_vx_A<<'\n' , file_vxres<<'\n', file_vy<<'\n', file_vy_A<<'\n' , file_vyres<<'\n')
            {
                file_vx    << xss(i, j) << ' ' << yss(i, j) << ' ' << vxss(i+nx*j);
                file_vy    << xss(i, j) << ' ' << yss(i, j) << ' ' << vyss(i+nx*j);
                file_vx_A  << xss(i, j) << ' ' << yss(i, j) << ' ' << vx_A(i+nx*j);
                file_vy_A  << xss(i, j) << ' ' << yss(i, j) << ' ' << vy_A(i+nx*j);
                file_vxres << xss(i, j) << ' ' << yss(i, j) << ' ' << vxss(i, j)-vx_A(i+nx*j);
                file_vyres << xss(i, j) << ' ' << yss(i, j) << ' ' << vyss(i, j)-vy_A(i+nx*j);
            }
        file_vx.close();
        file_vy.close();
        file_vx_A.close();
        file_vy_A.close();
        file_vxres.close();
        file_vyres.close();

        double res_1 = sqrt(sum(res_vx%res_vx)/sum(vx_A%vx_A));
        double res_2 = sqrt(sum(res_vy%res_vy)/sum(vy_A%vy_A));
        std::cout << "nx = " << nx << " ny = " << ny << " res_v1 = " << res_1 << " res_v2 = " << res_2 << '\n';

        file << nx << ' ' << ny << ' ' << res_1 << ' ' << res_2 << '\n';
    }
    file.close();
}