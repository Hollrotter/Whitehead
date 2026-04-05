#include "Membrane.hpp"
#include "interpolation.hpp"
#include <ginac/ginac.h>

using namespace GiNaC;

int main()
{
    symbol x1("x1"), x2("x2");

    size_t n = 21;
    size_t samples = 1;

    const double nu = 0.3;
    const double Et = 1000;
    const double D = Et/(1 - nu*nu);
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
    double a3 = 0.25;
    double a4 = 0.2;
    double a5 = 0.5;

    ex x = x1*(1 + a1*sin(pi*x2));
    ex y = x2*(1 + a2*sin(pi*x1));

    // Displacement
    // ex v_x = sin(x1)*cos(pi/8*y);// + 0.1*sin(pi*x2)*cos(pi/8*x2);
    // ex v_y = sin(x2)*cos(pi/8*x);// + 0.1*sin(pi*x1)*cos(pi/8*x1);
    ex v_x = a3*sin(pi/2*x)*cos(pi/4*y);
    ex v_y = a4*sin(pi/2*y)*cos(pi/4*x);
    // ex z = 0.1*(1 - pow(x1-1, 2))*(1 - pow(x2-1, 2));
    // ex z   = 0.2*(1 - pow(x, 2) - pow(y, 2));
    // ex z = 0.45*(1 - x1*x1)*(1 - x2*x2);
    // ex z = 0.1*(1 - (x1+1)*(x1+1)/4)*(1 - (x2+1)*(x2+1)/4);
    ex z = a5*cos(x)*cos(y);

    ex dxdx1 = diff(x, x1);
    ex dxdx2 = diff(x, x2);
    ex dydx1 = diff(y, x1);
    ex dydx2 = diff(y, x2);

    ex v_1 = dxdx1*v_x + dydx1*v_y;
    ex v_2 = dxdx2*v_x + dydx2*v_y;

    ex v_1_1 = diff(v_1, x1);
    ex v_1_2 = diff(v_1, x2);
    ex v_2_1 = diff(v_2, x1);
    ex v_2_2 = diff(v_2, x2);
    
    ex v_1_11 = diff(v_1_1, x1);
    ex v_1_12 = diff(v_1_1, x2);
    ex v_1_22 = diff(v_1_2, x2);
    ex v_2_11 = diff(v_2_1, x1);
    ex v_2_12 = diff(v_2_1, x2);
    ex v_2_22 = diff(v_2_2, x2);
    
    ex z_1  = diff(z,   x1);
    ex z_2  = diff(z,   x2);
    ex z_11 = diff(z_1, x1);
    ex z_12 = diff(z_1, x2);
    ex z_22 = diff(z_2, x2);

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
    
    ex v_1__11 = v_1_11 - 3*g111*v_1_1 - g211*(v_1_2 + 2*v_2_1)
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
    ex v_2__22 = v_2_22 - 3*g222*v_2_2 - g122*(v_2_1 + 2*v_1_2)
               + (2*g122*g112 + 2*g222*g122 - g122_2)*v_1
               + (2*g122*g212 + 2*g222*g222 - g222_2)*v_2;
    
    ex z__11 = z_11 - g111*z_1 - g211*z_2;
    ex z__12 = z_12 - g112*z_1 - g212*z_2;
    ex z__22 = z_22 - g122*z_1 - g222*z_2;
    
    ex g_11 =  v_1__1          + (z_1*z_1 + e11*v_1__1*v_1__1 + 2*e12*v_1__1*v_2__1                   + e22*v_2__1*v_2__1)/2;
    ex g_12 = (v_1__2 + v_2__1 +  z_1*z_2 + e11*v_1__1*v_1__2 +   e12*(v_1__1*v_2__2 + v_1__2*v_2__1) + e22*v_2__1*v_2__2)/2;
    ex g_22 =  v_2__2          + (z_2*z_2 + e11*v_1__2*v_1__2 + 2*e12*v_1__2*v_2__2                   + e22*v_2__2*v_2__2)/2;
    
    ex g_11_1 = diff(g_11, x1);
    ex g_11_2 = diff(g_11, x2);
    ex g_12_1 = diff(g_12, x1);
    ex g_12_2 = diff(g_12, x2);
    ex g_22_1 = diff(g_22, x1);
    ex g_22_2 = diff(g_22, x2);

    ex g_11__1 = g_11_1 - 2*g111*g_11 - 2*g211*g_12;
    ex g_11__2 = g_11_2 - 2*g112*g_11 - 2*g212*g_12;
    ex g_12__1 = g_12_1 -   g112*g_11 - (g111 + g212)*g_12 - g211*g_22;
    ex g_12__2 = g_12_2 -   g122*g_11 - (g112 + g222)*g_12 - g212*g_22;
    ex g_22__1 = g_22_1 - 2*g112*g_12 - 2*g212*g_22;
    ex g_22__2 = g_22_2 - 2*g122*g_12 - 2*g222*g_22;
    
    // ex a_11 = e_11 + 2*g_11;
    // ex a_12 = e_12 + 2*g_12;
    // ex a_22 = e_22 + 2*g_22;
    // ex a = a_11*a_22 - a_12*a_12;
    // ex ae = a/e;
    // ex sae = sqrt(ae);
    ex sae = 1 + e11*g_11 + 2*e12*g_12 + e22*g_22;

    ex H1111 = e11*e11;
    ex H1112 = e11*e12;
    ex H1122 = e12*e12 + nu/e;
    ex H1221 = (e11*e22 + e12*e12 - nu/e)/2;
    ex H1222 = e12*e22;
    ex H2222 = e22*e22;

    ex n11 = D*(H1111*g_11 + 2*H1112*g_12 + H1122*g_22);
    ex n12 = D*(H1112*g_11 + 2*H1221*g_12 + H1222*g_22);
    ex n22 = D*(H1122*g_11 + 2*H1222*g_12 + H2222*g_22);
    
    ex n11__1 = D*(H1111*g_11__1 + 2*H1112*g_12__1 + H1122*g_22__1);
    ex n12__1 = D*(H1112*g_11__1 + 2*H1221*g_12__1 + H1222*g_22__1);
    ex n12__2 = D*(H1112*g_11__2 + 2*H1221*g_12__2 + H1222*g_22__2);
    ex n22__2 = D*(H1122*g_11__2 + 2*H1222*g_12__2 + H2222*g_22__2);
    
    ex L111 = e11*v_1__11 + e12*v_2__11;
    ex L112 = e11*v_1__12 + e12*v_2__12;
    ex L122 = e11*v_1__22 + e12*v_2__22;
    ex L211 = e12*v_1__11 + e22*v_2__11;
    ex L212 = e12*v_1__12 + e22*v_2__12;
    ex L222 = e12*v_1__22 + e22*v_2__22;
    
    ex w   =-(e11*z_1*z_1 + 2*e12*z_1*z_2 + e22*z_2*z_2)/2;
    ex w_1 = z_1*(e11*v_1__1 + e12*v_2__1 - 1) + z_2*(e12*v_1__1 + e22*v_2__1);
    ex w_2 = z_1*(e11*v_1__2 + e12*v_2__2)     + z_2*(e12*v_1__2 + e22*v_2__2 - 1);

    ex w1 = e11*w_1 + e12*w_2;
    ex w2 = e12*w_1 + e22*w_2;

    ex wx = dxdx1*w1 + dxdx2*w2;
    ex wy = dydx1*w1 + dydx2*w2;

    ex s1 =-(n11__1 + n12__2 + L111*n11 + 2*L112*n12 + L122*n22)/sae;
    ex s2 =-(n12__1 + n22__2 + L211*n11 + 2*L212*n12 + L222*n22)/sae;
    ex s  =-(z__11*n11 + 2*z__12*n12 + z__22*n22)/sae;

    ex sx = dxdx1*s1 + dxdx2*s2;
    ex sy = dydx1*s1 + dydx2*s2;
    
    // ex nxx = dxdx1*dxdx1*n11 +               2*dxdx1*dxdx2*n12 + dxdx2*dxdx2*n22;
    // ex nxy = dxdx1*dydx1*n11 + (dxdx1*dydx2 + dxdx2*dydx1)*n12 + dxdx2*dydx2*n22;
    // ex nyy = dydx1*dydx1*n11 +               2*dydx1*dydx2*n12 + dydx2*dydx2*n22;
        
    ex det = dxdx1*dydx2 - dxdx2*dydx1;
    ex z_x =( dydx2*z_1 - dydx1*z_2)/det;
    ex z_y =(-dxdx2*z_1 + dxdx1*z_2)/det;

    ex det2 = det*det;
    ex v_x__x = (dydx2*dydx2*v_1__1 - dydx2*dydx1*(v_1__2 + v_2__1)           + dydx1*dydx1*v_2__2)/det2;
    ex v_x__y =-(dydx2*dxdx2*v_1__1 - dydx2*dxdx1*v_1__2 - dydx1*dxdx2*v_2__1 + dydx1*dxdx1*v_2__2)/det2;
    ex v_y__x =-(dxdx2*dydx2*v_1__1 - dxdx2*dydx1*v_1__2 - dxdx1*dydx2*v_2__1 + dxdx1*dydx1*v_2__2)/det2;
    ex v_y__y = (dxdx2*dxdx2*v_1__1 - dxdx2*dxdx1*(v_1__2 + v_2__1)           + dxdx1*dxdx1*v_2__2)/det2;

    Membrane m;
    m.iterations(100);
    m.youngsModulus(Et);
    m.poissonsRatio(nu);
    std::ofstream file("plot/Data/MembraneNonlinearMMS");

    for (size_t k = 0; k < samples; k++, n+=2)
    {
        size_t n2 = n*n;
        arma::vec x1s = Chebyshev::gaussLobatto(n);
        arma::vec x2s = Chebyshev::gaussLobatto(n);
        
        arma::mat  xs(n, n), ys(n, n);
        arma::vec sxs(n2), sys(n2), ss(n2);
        arma::vec wxs(n2), wys(n2), ws(n2);
        arma::vec vx_A(n2), vy_A(n2), z_A(n2), zx_A(n2), zy_A(n2);
        // arma::vec nxx_A(n2), nxy_A(n2), nyy_A(n2);
        arma::vec v_x__x_A(n2), v_x__y_A(n2), v_y__x_A(n2), v_y__y_A(n2);

        exmap mp;
        std::cout << "j = " << std::flush;
        for (size_t j = 0; j < n; j++)
        {
            std::cout << j << ' ' << std::flush;
            mp[x2] = x2s(j);
            for (size_t i = 0; i < n; i++)
            {
                size_t q = i+j*n;
                mp[x1] = x1s(i);
                xs(i, j)    = ex_to<numeric>(evalf(  x.subs(mp))).to_double();
                ys(i, j)    = ex_to<numeric>(evalf(  y.subs(mp))).to_double();
                sxs(q)      = ex_to<numeric>(evalf( sx.subs(mp))).to_double();
                sys(q)      = ex_to<numeric>(evalf( sy.subs(mp))).to_double();
                ss(q)       = ex_to<numeric>(evalf(  s.subs(mp))).to_double();
                wxs(q)      = ex_to<numeric>(evalf( wx.subs(mp))).to_double();
                wys(q)      = ex_to<numeric>(evalf( wy.subs(mp))).to_double();
                ws(q)       = ex_to<numeric>(evalf(  w.subs(mp))).to_double();
                vx_A(q)     = ex_to<numeric>(evalf(v_x.subs(mp))).to_double();
                vy_A(q)     = ex_to<numeric>(evalf(v_y.subs(mp))).to_double();
                // nxx_A(q)    = ex_to<numeric>(evalf(nxx.subs(mp))).to_double();
                // nxy_A(q)    = ex_to<numeric>(evalf(nxy.subs(mp))).to_double();
                // nyy_A(q)    = ex_to<numeric>(evalf(nyy.subs(mp))).to_double();
                v_x__x_A(q) = ex_to<numeric>(evalf(v_x__x.subs(mp))).to_double();
                v_x__y_A(q) = ex_to<numeric>(evalf(v_x__y.subs(mp))).to_double();
                v_y__x_A(q) = ex_to<numeric>(evalf(v_y__x.subs(mp))).to_double();
                v_y__y_A(q) = ex_to<numeric>(evalf(v_y__y.subs(mp))).to_double();
                z_A(q)      = ex_to<numeric>(evalf(  z.subs(mp))).to_double();
                zx_A(q)     = ex_to<numeric>(evalf(z_x.subs(mp))).to_double();
                zy_A(q)     = ex_to<numeric>(evalf(z_y.subs(mp))).to_double();
            }
        }
        std::cout << '\n';
        Lagrange::CurveInterpolant chi1(xs.col(0),       ys.col(0));
        Lagrange::CurveInterpolant chi2(xs.row(n-1).t(), ys.row(n-1).t());
        Lagrange::CurveInterpolant chi3(xs.col(n-1),     ys.col(n-1));
        Lagrange::CurveInterpolant chi4(xs.row(0).t(),   ys.row(0).t());
        
        auto [xss, yss] = Lagrange::TransfiniteQuadMap({&chi1, &chi2, &chi3, &chi4});
        arma::mat T = ChebyshevInterpolation(x1s, x2s, xs, ys, xss, yss);
        vx_A      = T * vx_A;
        vy_A      = T * vy_A;
        z_A       = T * z_A;
        zx_A      = T * zx_A;
        zy_A      = T * zy_A;
        // nxx_A     = T * nxx_A;
        // nxy_A     = T * nxy_A;
        // nyy_A     = T * nyy_A;
        v_x__x_A  = T * v_x__x_A;
        v_x__y_A  = T * v_x__y_A;
        v_y__x_A  = T * v_y__x_A;
        v_y__y_A  = T * v_y__y_A;

        arma::vec sxss = T * sxs;
        arma::vec syss = T * sys;
        arma::vec sss  = T * ss;

        arma::vec wxss = T * wxs;
        arma::vec wyss = T * wys;
        arma::vec wss  = T * ws;

        m({&chi1, &chi2, &chi3, &chi4});

        arma::mat e_11ss = m.e_11();
        arma::mat e_12ss = m.e_12();
        arma::mat e_22ss = m.e_22();
        arma::mat e11ss = m.e11();
        arma::mat e12ss = m.e12();
        arma::mat e22ss = m.e22();

        arma::mat dxdx1ss = m.J11();
        arma::mat dxdx2ss = m.J12();
        arma::mat dydx1ss = m.J21();
        arma::mat dydx2ss = m.J22();

        arma::mat det_ss = dxdx1ss%dydx2ss - dydx1ss%dxdx2ss;
        arma::mat det2_ss = det_ss%det_ss;

        arma::mat v1_A  = dxdx1ss%reshape(vx_A, n, n) + dydx1ss%reshape(vy_A, n, n);
        arma::mat v2_A  = dxdx2ss%reshape(vx_A, n, n) + dydx2ss%reshape(vy_A, n, n);

        // arma::mat n11_A = (dydx2ss%dydx2ss%reshape(nxx_A, n, n) -                   2*dydx2ss%dxdx2ss%reshape(nxy_A, n, n) + dxdx2ss%dxdx2ss%reshape(nyy_A, n, n))/det2_ss;
        // arma::mat n12_A =-(dydx2ss%dydx1ss%reshape(nxx_A, n, n) - (dydx2ss%dxdx1ss + dydx1ss%dxdx2ss)%reshape(nxy_A, n, n) + dxdx2ss%dxdx1ss%reshape(nyy_A, n, n))/det2_ss;
        // arma::mat n22_A = (dydx1ss%dydx1ss%reshape(nxx_A, n, n) -                   2*dydx1ss%dxdx1ss%reshape(nxy_A, n, n) + dxdx1ss%dxdx1ss%reshape(nyy_A, n, n))/det2_ss;

        arma::mat w1_A = ( dydx2ss%reshape(wxss, n, n) - dxdx2ss%reshape(wyss, n, n))/det_ss;
        arma::mat w2_A = (-dydx1ss%reshape(wxss, n, n) + dxdx1ss%reshape(wyss, n, n))/det_ss;

        arma::mat s1_A = ( dydx2ss%reshape(sxss, n, n) - dxdx2ss%reshape(syss, n, n))/det_ss;
        arma::mat s2_A = (-dydx1ss%reshape(sxss, n, n) + dxdx1ss%reshape(syss, n, n))/det_ss;

        arma::mat z_1_A = dxdx1ss%reshape(zx_A, n, n) + dydx1ss%reshape(zy_A, n, n);
        arma::mat z_2_A = dxdx2ss%reshape(zx_A, n, n) + dydx2ss%reshape(zy_A, n, n);

        arma::mat v_1__1_A = dxdx1ss%dxdx1ss%reshape(v_x__x_A, n, n) + dxdx1ss%dydx1ss%reshape(v_x__y_A + v_y__x_A, n, n)                                + dydx1ss%dydx1ss%reshape(v_y__y_A, n, n);
        arma::mat v_1__2_A = dxdx1ss%dxdx2ss%reshape(v_x__x_A, n, n) + dxdx1ss%dydx2ss%reshape(v_x__y_A, n, n) + dydx1ss%dxdx2ss%reshape(v_y__x_A, n, n) + dydx1ss%dydx2ss%reshape(v_y__y_A, n, n);
        arma::mat v_2__1_A = dxdx2ss%dxdx1ss%reshape(v_x__x_A, n, n) + dxdx2ss%dydx1ss%reshape(v_x__y_A, n, n) + dydx2ss%dxdx1ss%reshape(v_y__x_A, n, n) + dydx2ss%dydx1ss%reshape(v_y__y_A, n, n);
        arma::mat v_2__2_A = dxdx2ss%dxdx2ss%reshape(v_x__x_A, n, n) + dxdx2ss%dydx2ss%reshape(v_x__y_A + v_y__x_A, n, n)                                + dydx2ss%dydx2ss%reshape(v_y__y_A, n, n);

        // Boundary Conditions

        // West and East
        arma::mat h_11s = sqrt(1/e11ss);
        arma::mat h_12s = e_12ss/sqrt(e_22ss);
        arma::mat h_22s = sqrt(e_22ss);
        
        arma::mat h_1s1 = sqrt(e11ss);
        arma::mat h_1s2 = e12ss/sqrt(e11ss);
        arma::mat h_2s2 = 1/sqrt(e_22ss);
        
        arma::mat v_1s_WE = h_1s1%v1_A + h_1s2%v2_A;
        arma::mat v_2s_WE = h_2s2%v2_A;
        arma::mat z_1s_WE = h_1s1%z_1_A + h_1s2%z_2_A;
        // arma::mat v_1s__1s = h_1s1*h_1s1*v_1__1 + h_1s1*h_1s2*(v_1__2 + v_2__1) + h_1s2*h_1s2*v_2__2;
        // arma::mat v_2s__1s = h_1s1*h_2s2*v_2__1 + h_2s2*h_1s2*v_2__2;
        
        // arma::mat n_1s1s_WE = h_11s%h_11s%n11_A;
        // arma::mat n_1s2s_WE = h_11s%h_12s%n11_A + h_11s%h_22s%n12_A;

        // North and South
        h_11s    = sqrt(e_11ss);
        arma::mat h_21s = e_12ss/sqrt(e_11ss);
        h_22s    = sqrt(1/e22ss);
        
        h_1s1    = 1/sqrt(e_11ss);
        arma::mat h_2s1 = e12ss/sqrt(e22ss);
        h_2s2    = sqrt(e22ss);
        
        arma::mat v_1s_NS = h_1s1%v1_A;
        arma::mat v_2s_NS = h_2s1%v1_A + h_2s2%v2_A;
        arma::mat z_2s_NS = h_2s1%z_1_A + h_2s2%z_2_A;
        // arma::mat v_1s__2s = h_1s1*h_2s1*v_1__1 + h_1s1*h_2s2*v_1__2;
        // arma::mat v_2s__2s = h_2s1*h_2s1*v_1__1 + h_2s1*h_2s2*(v_1__2 + v_2__1) + h_2s2*h_2s2*v_2__2;
        
        // arma::mat n_2s2s_NS = h_22s%h_22s%n22_A;
        // arma::mat n_2s1s_NS = h_11s%h_22s%n12_A + h_21s%h_22s%n22_A;

        arma::vec v1_N(n), v1_S(n), v1_W(n), v1_E(n);
        arma::vec v2_N(n), v2_S(n), v2_W(n), v2_E(n);
        arma::vec  z_N(n),  z_S(n),  z_W(n),  z_E(n);
        // arma::vec n11_W(ny), n11_E(ny);
        // arma::vec n12_W(ny), n12_E(ny);
        // arma::vec n21_N(nx), n21_S(nx);
        // arma::vec n22_N(nx), n22_S(nx);

        for (size_t i = 0; i < n; i++)
        {
            v1_S(i)  = v_1s_NS(i,   0);
            v2_S(i)  = v_2s_NS(i,   0);
            v1_N(i)  = v_1s_NS(i, n-1);
            v2_N(i)  = v_2s_NS(i, n-1);
            z_S(i)   = r1_S*z_A(i)           + r2_S*z_2s_NS(i,   0);
            z_N(i)   = r1_N*z_A(i+n*(n-1)) + r2_N*z_2s_NS(i, n-1);
            // n22_S(i) = n_2s2s_NS(i, 0);
            // n21_S(i) = n_2s1s_NS(i, 0);
            // n22_N(i) = n_2s2s_NS(i, ny-1);
            // n21_N(i) = n_2s1s_NS(i, ny-1);
        }

        for (size_t j = 0; j < n; j++)
        {
            v1_W(j)  = v_1s_WE(  0, j);
            v2_W(j)  = v_2s_WE(  0, j);
            v1_E(j)  = v_1s_WE(n-1, j);
            v2_E(j)  = v_2s_WE(n-1, j);
            z_W(j)   = r1_W*z_A(n*j)     + r2_W*z_1s_WE(  0, j);
            z_E(j)   = r1_E*z_A(n*j+n-1) + r2_E*z_1s_WE(n-1, j);
            // n11_W(j) = n_1s1s_WE(0, j);
            // n12_W(j) = n_1s2s_WE(0, j);
            // n11_E(j) = n_1s1s_WE(nx-1, j);
            // n12_E(j) = n_1s2s_WE(nx-1, j);
        }

        m.boundary(Field::v1, Direction::N, BC::Dirichlet, v1_N);
        m.boundary(Field::v1, Direction::S, BC::Dirichlet, v1_S);
        m.boundary(Field::v1, Direction::W, BC::Dirichlet, v1_W);
        m.boundary(Field::v1, Direction::E, BC::Dirichlet, v1_E);

        // m.boundary(Field::n12, Direction::N, BC::Dirichlet, n21_N);
        // m.boundary(Field::n12, Direction::S, BC::Dirichlet, n21_S);
        // m.boundary(Field::n11, Direction::W, BC::Dirichlet, n11_W);
        // m.boundary(Field::n11, Direction::E, BC::Dirichlet, n11_E);

        m.boundary(Field::v2, Direction::N, BC::Dirichlet, v2_N);
        m.boundary(Field::v2, Direction::S, BC::Dirichlet, v2_S);
        m.boundary(Field::v2, Direction::W, BC::Dirichlet, v2_W);
        m.boundary(Field::v2, Direction::E, BC::Dirichlet, v2_E);

        // m.boundary(Field::n22, Direction::N, BC::Dirichlet, n22_N);
        // m.boundary(Field::n22, Direction::S, BC::Dirichlet, n22_S);
        // m.boundary(Field::n12, Direction::W, BC::Dirichlet, n12_W);
        // m.boundary(Field::n12, Direction::E, BC::Dirichlet, n12_E);

        m.boundary(Field::z,  Direction::N, BC::Robin, r1_N, r2_N, z_N);
        m.boundary(Field::z,  Direction::S, BC::Robin, r1_S, r2_S, z_S);
        m.boundary(Field::z,  Direction::W, BC::Robin, r1_W, r2_W, z_W);
        m.boundary(Field::z,  Direction::E, BC::Robin, r1_E, r2_E, z_E);

        arma::mat pss(n, n);
        arma::vec p1_A(n2), p2_A(n2);
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < n; i++)
            {
                arma::mat M = {{e11ss(i, j)*v_1__1_A(i, j) + e12ss(i, j)*v_2__1_A(i, j) + 1, e11ss(i, j)*v_1__2_A(i, j) + e12ss(i, j)*v_2__2_A(i, j),     w1_A(i, j)},
                               {e12ss(i, j)*v_1__1_A(i, j) + e22ss(i, j)*v_2__1_A(i, j)    , e12ss(i, j)*v_1__2_A(i, j) + e22ss(i, j)*v_2__2_A(i, j) + 1, w2_A(i, j)},
                               {z_1_A(i, j), z_2_A(i, j), 1 + wss(i+j*n)}};
                arma::vec rhs = {s1_A(i+j*n), s2_A(i+j*n), sss(i+j*n)};
                arma::vec ppp = solve(M, rhs);
                p1_A(i+j*n) = ppp(0);
                p2_A(i+j*n) = ppp(1);
                pss(i, j)   = ppp(2);
            }
        m.inPlaneTensorial(p1_A, p2_A);
        m.load(pss);
        printf("Starting solution...\n");
        if (k == 0)
            m.planeStrain();
        m.nonlinear();

        auto [vxss, vyss] = m.kartesianDisplacements();

        arma::vec z_num = vectorise(arma::mat(m(Field::z)));

        m.output("plot/Data/z",  Field::z);
        m.output("plot/Data/v1", Field::v1);
        m.output("plot/Data/v2", Field::v2);
        std::ofstream file_z("plot/Data/z_A");
        std::ofstream file_vx("plot/Data/vx");
        std::ofstream file_vy("plot/Data/vy");
        std::ofstream file_vx_A("plot/Data/vx_A");
        std::ofstream file_vy_A("plot/Data/vy_A");
        std::ofstream file_vxres("plot/Data/vx_res");
        std::ofstream file_vyres("plot/Data/vy_res");
        std::ofstream file_zres("plot/Data/z_res");
        std::ofstream file_xys ("plot/Data/xys");
        std::ofstream file_xyss("plot/Data/xyss");
        for (size_t j = 0; j < n; j++, file_z<<'\n', file_zres<<'\n', file_vx<<'\n', file_vx_A<<'\n' , file_vxres<<'\n', file_vy<<'\n', file_vy_A<<'\n' , file_vyres<<'\n', file_xys<<'\n', file_xyss<<'\n')
            for (size_t i = 0; i < n; i++, file_z<<'\n', file_zres<<'\n', file_vx<<'\n', file_vx_A<<'\n' , file_vxres<<'\n', file_vy<<'\n', file_vy_A<<'\n' , file_vyres<<'\n', file_xys<<'\n', file_xyss<<'\n')
            {
                file_z     << xss(i, j) << ' ' << yss(i, j) << ' ' << z_A(i+n*j);
                file_zres  << xss(i, j) << ' ' << yss(i, j) << ' ' << z_num(i+n*j)-z_A(i+n*j);
                file_vx    << xss(i, j) << ' ' << yss(i, j) << ' ' << vxss(i+n*j);
                file_vy    << xss(i, j) << ' ' << yss(i, j) << ' ' << vyss(i+n*j);
                file_vx_A  << xss(i, j) << ' ' << yss(i, j) << ' ' << vx_A(i+n*j);
                file_vy_A  << xss(i, j) << ' ' << yss(i, j) << ' ' << vy_A(i+n*j);
                file_vxres << xss(i, j) << ' ' << yss(i, j) << ' ' << vxss(i, j)-vx_A(i+n*j);
                file_vyres << xss(i, j) << ' ' << yss(i, j) << ' ' << vyss(i, j)-vy_A(i+n*j);
                file_xys  << xs(i, j)  << ' ' << ys(i, j)  << " 0";
                file_xyss << xss(i, j) << ' ' << yss(i, j) << " 0";
            }
        file_z.close();
        file_vx.close();
        file_vy.close();
        file_vx_A.close();
        file_vy_A.close();
        file_vxres.close();
        file_vyres.close();
        file_zres.close();
        file_xys.close();
        file_xyss.close();

        arma::vec res_vx = vectorise(vxss) - vx_A;
        arma::vec res_vy = vectorise(vyss) - vy_A;
        arma::vec res_z  = z_num - z_A;
        double res_1 = sqrt(sum(res_vx % res_vx)/sum(vx_A % vx_A));
        double res_2 = sqrt(sum(res_vy % res_vy)/sum(vy_A % vy_A));
        double res_3 = sqrt(sum(res_z  % res_z) /sum( z_A %  z_A));
        std::cout << "n = " << n << " res_v1 = " << res_1 << " res_v2 = " << res_2 << " res_z = " << res_3 << '\n';

        file << n << ' ' << res_1 << ' ' << res_2 << ' ' << res_3 << '\n';
    }
    file.close();
}