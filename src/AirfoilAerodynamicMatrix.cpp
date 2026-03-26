#include "Airfoil.hpp"

/**
 * @brief 
 * 
 */
void Airfoil::aerodynamicMatrix()
{
    switch (analysis)
    {
        case Analysis::linear:
        {
            for (size_t i = 0; i < nx-1; i++)
            {
                A(i, 0) = log((1 - xi(i))/(1 + xi(i)));
                A(i, 1) = 2 + xi(i)*A(i, 0);
                for (size_t j = 2; j < nx; j++)
                    A(i, j) = 2*(Chebyshev::integral(j-1) + xi(i)*A(i, j-1)) - A(i, j-2);
                b.row(i) =-arma::datum::tau*alpha.t();
            }
            A.row(nx-1).ones();
            break;
        }
        case Analysis::nonlinear:
        {
            size_t mx = 2;
            std::vector<fastgl::QuadPair> gl_n(nx);
            for (size_t n = 0; n < nx; n++)
                gl_n[n] = fastgl::GLPair(nx, n+1);
            std::vector<fastgl::QuadPair> gl_m(mx);
            for (size_t m = 0; m < mx; m++)
                gl_m[m] = fastgl::GLPair(mx, m+1);
            nC.zeros(2, nx);
            auto [dxdxi, dzdxi] = chi->derivative(xi);
            arma::mat D = Lagrange::derivativeMatrix(xi);
            arma::vec d2xdxi2 = D*dxdxi;
            arma::vec d2zdxi2 = D*dzdxi;
            arma::mat A_t = join_horiz(dxdxi, dzdxi);
            arma::mat B_t = join_horiz(d2xdxi2, d2zdxi2)/2;
            arma::vec C_t = sum(A_t%B_t, 1);
            for (size_t i = 0; i < nx-1; i++)
            {
                double xi_min = i == 0    ? -1 : xi(i-1);
                double xi_max = i == nx-1 ?  1 : xi(i+1);
                double J = sqrt(pow(dxdxi(i), 2) + pow(dzdxi(i), 2));
                nC.col(i) = arma::vec{-dzdxi(i), dxdxi(i)}/J;
                for (size_t j = 0; j < nx; j++)
                {
                    double gamma_0 = boost::math::chebyshev_t(j, xi(i));
                    double gamma_1 = boost::math::chebyshev_t_prime(j, xi(i));
                    double A_abs = norm(A_t.row(i));
                    double A4 = pow(A_abs, 4);
                    arma::vec::fixed<2> c1 = A_t.row(i).t()/A_abs;
                    arma::vec::fixed<2> c2 = B_t.row(i).t()/A_abs - 2*C_t(i)/A4;
                    arma::vec::fixed<2> c3 =-2*C_t(i)/A4*B_t.row(i).t();
                    arma::vec::fixed<2> F_1 = gamma_0 * arma::vec::fixed<2>{c1(1),-c1(0)};
                    arma::vec::fixed<2> F0  = gamma_0 * arma::vec::fixed<2>{c2(1),-c2(0)} + gamma_1 * arma::vec::fixed<2>{c1(1),-c1(0)};
                    arma::vec::fixed<2> F1  = gamma_0 * arma::vec::fixed<2>{c3(1),-c3(0)} + gamma_1 * arma::vec::fixed<2>{c2(1),-c2(0)};
                    arma::vec::fixed<2> F2  = gamma_1 * arma::vec::fixed<2>{c3(1),-c3(0)};
                    arma::vec::fixed<2> I_r(arma::fill::zeros);
                    for (size_t m = 0; m < mx; m++)
                    {
                        double rho = (xi_max - xi_min)/2*(1 - gl_m[m].x());
                        I_r += gl_m[m].weight * (F0 + F1*rho + F2*pow(rho, 2));
                    }
                    arma::vec::fixed<2> q = I_r * (xi_max - xi_min)/2;
                    q += F_1*log((1-xi(i))/(1+xi(i))) * (xi_max - xi_min)/2;
                    if (i > 0)
                    {
                        arma::vec::fixed<2> I;
                        for (size_t n = 0; n < nx; n++)
                        {
                            auto [x_n, z_n] = chi->evaluate((1+gl_n[n].x())/2*(1-xi_max));
                            arma::vec::fixed<2> r ={z(i)-z_n, x_n-x(i)};
                            double r2 = dot(r, r);
                            double T = boost::math::chebyshev_t(j, (1+gl_n[n].x())/2*(1-xi_max));
                            I += gl_n[n].weight * T * r/r2;
                        }
                        q += I * (1-xi_max)/2;
                    }
                    if (i < nx-1)
                    {
                        arma::vec::fixed<2> I;
                        for (size_t n = 0; n < nx; n++)
                        {
                            auto [x_n, z_n] = chi->evaluate((1+gl_n[n].x())/2*(xi_min-1));
                            arma::vec::fixed<2> r = {z(i)-z_n, x_n-x(i)};
                            double r2 = dot(r, r);
                            double T = boost::math::chebyshev_t(j, (1+gl_n[n].x())/2*(xi_min-1));
                            I += gl_n[n].weight * T * r/r2;
                        }
                        q += I * (xi_min-1)/2;
                    }
                    A(i, j) = dot(nC.col(i), q);
                }
            }
            A.row(nx-1).ones();
            break;
        }
    }
}