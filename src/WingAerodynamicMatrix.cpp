#include "Wing.hpp"

arma::vec Wing::externalContour(double x_min, double x_max, double y_min, double y_max, double eta_1, double eta_2, arma::vec theta)
{
    arma::vec rho(theta.size());
    for (size_t n = 0; n < theta.size(); n++)
        if (-arma::datum::pi < theta(n) && theta(n) <=-arma::datum::pi/2 || arma::datum::pi < theta(n) && theta(n) <= 3*arma::datum::pi/2)
            rho(n) = std::min(fabs((x_min-eta_1)/cos(theta(n))), fabs((y_min-eta_2)/sin(theta(n))));
        else if (-arma::datum::pi/2 < theta(n) && theta(n) <= 0 || 3*arma::datum::pi/2 < theta(n) && theta(n) <= 2*arma::datum::pi)
            rho(n) = std::min(fabs((x_max-eta_1)/cos(theta(n))), fabs((y_min-eta_2)/sin(theta(n))));
        else if (0 < theta(n) && theta(n) <= arma::datum::pi/2)
            rho(n) = std::min(fabs((x_max-eta_1)/cos(theta(n))), fabs((y_max-eta_2)/sin(theta(n))));
        else
            rho(n) = std::min(fabs((x_min-eta_1)/cos(theta(n))), fabs((y_max-eta_2)/sin(theta(n))));
    return rho;
}

/**
 * @brief 
 * 
 */
void Wing::aerodynamicMatrix()
{
    size_t n_theta = 20;
    arma::vec theta(n_theta);
    std::vector<fastgl::QuadPair> gl_theta(n_theta);
    for (size_t n = 0; n < n_theta; n++)
    {
        gl_theta[n] = fastgl::GLPair(n_theta, n+1);
        theta(n) =-arma::datum::pi * gl_theta[n].x();
    }
    arma::vec cT = cos(theta);
    arma::vec sT = sin(theta);
    arma::vec scT = sin(theta)%cos(theta);
    arma::vec c2T = pow(cos(theta), 2)/2;
    arma::vec s2T = pow(sin(theta), 2)/2;
    switch (analysis)
    {
        case Analysis::linear:
        {
            size_t m_rho = 1;
            std::vector<fastgl::QuadPair> gl_rho(m_rho);
            for (size_t m = 0; m < m_rho; m++)
                gl_rho[m] = fastgl::GLPair(m_rho, m+1);
            auto [xC, yC] = Lagrange::TransfiniteQuadMap(xi_1, xi_2, chi);
            auto [dxdxi_1, dxdxi_2, dydxi_1, dydxi_2] = Lagrange::TransfiniteQuadMetrics(xi_1, xi_2, chi);
            arma::mat detJ = dxdxi_1%dydxi_2 - dxdxi_2%dydxi_1;
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double eta_2 = xi_2(j);
                double dxi_2_m = eta_2 - xi_2(j-1);
                double dxi_2_p = xi_2(j+1) - eta_2;
                double DY = dxi_2_m*dxi_2_p*(dxi_2_m + dxi_2_p);

                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double eta_1 = xi_1(i);
                    double dxi_1_m = eta_1 - xi_1(i-1);
                    double dxi_1_p = xi_1(i+1) - eta_1;
                    double DX = dxi_1_m*dxi_1_p*(dxi_1_m + dxi_1_p);

                    arma::vec::fixed<2> dXdxi_1 = {dxdxi_1(i, j), dydxi_1(i, j)};
                    arma::vec::fixed<2> dXdxi_2 = {dxdxi_2(i, j), dydxi_2(i, j)};

                    double d2xdxi_12     = (-pow(dxi_1_p, 2)*dxdxi_1(i-1, j) + (pow(dxi_1_p, 2) - pow(dxi_1_m, 2))*dxdxi_1(i, j) + pow(dxi_1_m, 2)*dxdxi_1(i+1, j))/DX;
                    double d2xdxi_1dxi_2 = (-pow(dxi_2_p, 2)*dxdxi_1(i, j-1) + (pow(dxi_2_p, 2) - pow(dxi_2_m, 2))*dxdxi_1(i, j) + pow(dxi_2_m, 2)*dxdxi_1(i, j+1))/DY;
                    double d2xdxi_22     = (-pow(dxi_2_p, 2)*dxdxi_2(i, j-1) + (pow(dxi_2_p, 2) - pow(dxi_2_m, 2))*dxdxi_2(i, j) + pow(dxi_2_m, 2)*dxdxi_2(i, j+1))/DY;
                    double d2ydxi_12     = (-pow(dxi_1_p, 2)*dydxi_1(i-1, j) + (pow(dxi_1_p, 2) - pow(dxi_1_m, 2))*dydxi_1(i, j) + pow(dxi_1_m, 2)*dydxi_1(i+1, j))/DX;
                    double d2ydxi_1dxi_2 = (-pow(dxi_1_p, 2)*dydxi_2(i-1, j) + (pow(dxi_1_p, 2) - pow(dxi_1_m, 2))*dydxi_2(i, j) + pow(dxi_1_m, 2)*dydxi_2(i+1, j))/DX;
                    double d2ydxi_22     = (-pow(dxi_2_p, 2)*dydxi_2(i, j-1) + (pow(dxi_2_p, 2) - pow(dxi_2_m, 2))*dydxi_2(i, j) + pow(dxi_2_m, 2)*dydxi_2(i, j+1))/DY;

                    arma::vec::fixed<2> d2Xdxi_12     = {d2xdxi_12,     d2ydxi_12};
                    arma::vec::fixed<2> d2Xdxi_1dxi_2 = {d2xdxi_1dxi_2, d2ydxi_1dxi_2};
                    arma::vec::fixed<2> d2Xdxi_22     = {d2xdxi_22,     d2ydxi_22};

                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, eta_1, eta_2, theta);

                    size_t k1 = i+j*nx;
                    b.row(k1) =-4*arma::datum::pi*alpha; // Must be corrected later!

                    double J0 = detJ(i, j);
                    double dJdxi_1 = (-pow(dxi_1_p, 2)*detJ(i-1, j) + (pow(dxi_1_p, 2) - pow(dxi_1_m, 2))*detJ(i, j) + pow(dxi_1_m, 2)*detJ(i+1, j))/DX;
                    double dJdxi_2 = (-pow(dxi_2_p, 2)*detJ(i, j-1) + (pow(dxi_2_p, 2) - pow(dxi_2_m, 2))*detJ(i, j) + pow(dxi_2_m, 2)*detJ(i, j+1))/DY;
                    arma::vec J1 = dJdxi_1*cT + dJdxi_2*sT;
                    
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T_1 = boost::math::chebyshev_t_prime(p, eta_1) * boost::math::chebyshev_t(q, eta_2);
                            double T_2 = boost::math::chebyshev_t(p, eta_1)       * boost::math::chebyshev_t_prime(p, eta_2);
                            double mu0 = boost::math::chebyshev_t(p, eta_1) * boost::math::chebyshev_t(q, eta_2);
                            double I0  = 0;
                            double I_1 = 0;
                            double I_2 = 0;
                            arma::vec mu1 = T_1*cT + T_2*sT;
                            for (size_t n = 0; n < n_theta; n++)
                            {
                                arma::vec::fixed<2> A_t = dXdxi_1*cT(n) + dXdxi_2*sT(n);
                                arma::vec::fixed<2> B_t = d2Xdxi_12*c2T(n) + d2Xdxi_1dxi_2*scT(n) + d2Xdxi_22*s2T(n);
                                double C_t = dot(A_t, B_t);
                                double A_abs = norm(A_t);
                                double A2 = pow(A_abs, 2);
                                double A3 = pow(A_abs, 3);
                                double A4 = pow(A2, 2);
                                double A5 = A2*A3;
                                double F_2 = J0*mu0/A3;
                                double F_1 = ((J1(n) - 3/A2*C_t*J0)*mu0 + J0*mu1(n))/A3;
                                double F0  =-3/A5*(C_t*J1(n)*mu0 + (C_t*J0 - J1(n)/3*A2)*mu1(n));
                                double F1  =-3/A5*J1(n)*C_t*mu1(n);

                                double I_r = 0;
                                for (size_t m = 0; m < m_rho; m++)
                                    I_r += gl_rho[m].weight * (F0 + F1*rho_tilde(n)*(1 - gl_rho[m].x())/2);
                                I_r *= rho_tilde(n)/2;
                                I0  += gl_theta[n].weight * I_r;
                                I_1 += gl_theta[n].weight * F_1*log(fabs(rho_tilde(n)*A_abs));
                                I_2 -= gl_theta[n].weight * F_2*(1/rho_tilde(n) + C_t/A2);
                            }
                            A(k1, p+q*nx) = arma::datum::pi*(I0 + I_1 + I_2);
                        }
                }
            }
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    size_t k1 = i + j*nx;
                    if (j > 1)
                    {
                        std::vector<fastgl::QuadPair> gl_x(nx), gl_y(j+2);
                        arma::vec x1_gl(nx), x2_gl(j+2);
                        for (size_t ii = 0; ii < nx; ii++)
                        {
                            gl_x[ii] = fastgl::GLPair(nx, ii+1);
                            x1_gl(ii) =-gl_x[ii].x();
                        }
                        for (size_t jj = 0; jj < j+2; jj++)
                        {
                            gl_y[jj] = fastgl::GLPair(j+2, jj+1);
                            x2_gl(jj) = (1 - gl_y[jj].x()) * (xi_2(j-1) + 1)/2 - 1;
                        }
                        auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                        auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                        arma::mat J_gl = dx_gldx1%dy_gldx2 - dx_gldx2%dy_gldx1;
                        for (size_t ii = 0; ii < nx; ii++)
                            for (size_t jj = 0; jj < j+2; jj++)
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        A(k1, p+q*nx) += gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                      / pow(sqrt(pow(x_gl(ii, jj) - xC(i, j), 2) + pow(y_gl(ii, jj) - yC(i, j), 2)), 3)
                                                      * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                      * (xi_2(j-1) + 1)/2;
                    }
                    if (j < ny-1)
                    {
                        std::vector<fastgl::QuadPair> gl_x(nx), gl_y(ny-j+2);
                        arma::vec x1_gl(nx), x2_gl(ny-j+2);
                        for (size_t ii = 0; ii < nx; ii++)
                        {
                            gl_x[ii] = fastgl::GLPair(nx, ii+1);
                            x1_gl(ii) =-gl_x[ii].x();
                        }
                        for (size_t jj = 0; jj < ny-j+2; jj++)
                        {
                            gl_y[jj] = fastgl::GLPair(ny-j+2, jj+1);
                            x2_gl(jj) = (1 - gl_y[jj].x()) * (1 - xi_2(j+1))/2 + xi_2(j+1);
                        }
                        auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                        auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                        arma::mat J_gl = dx_gldx1%dy_gldx2 - dx_gldx2%dy_gldx1;
                        for (size_t ii = 0; ii < nx; ii++)
                            for (size_t jj = 0; jj < ny-j+2; jj++)
                            {
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        A(k1, p+q*nx) += gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                      / pow(sqrt(pow(x_gl(ii, jj) - xC(i, j), 2) + pow(y_gl(ii, jj) - yC(i, j), 2)), 3)
                                                      * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                      * (1 - xi_2(j+1))/2;
                            }
                    }
                    if (i > 1)
                    {
                        std::vector<fastgl::QuadPair> gl_x(i+2), gl_y(4);
                        arma::vec x1_gl(i+2), x2_gl(4);
                        for (size_t ii = 0; ii < i+2; ii++)
                        {
                            gl_x[ii] = fastgl::GLPair(i+2, ii+1);
                            x1_gl(ii) = (1 - gl_x[ii].x()) * (xi_1(i-1) + 1)/2 - 1;
                        }
                        for (size_t jj = 0; jj < 4; jj++)
                        {
                            gl_y[jj] = fastgl::GLPair(4, jj+1);
                            x2_gl(jj) = (1 - gl_y[jj].x()) * (xi_2(j+1) - xi_2(j-1))/2 + xi_2(j-1);
                        }
                        auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                        auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                        arma::mat J_gl = dx_gldx1%dy_gldx2 - dx_gldx2%dy_gldx1;
                        for (size_t ii = 0; ii < i+2; ii++)
                            for (size_t jj = 0; jj < 4; jj++)
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        A(k1, p+q*nx) += gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                      / pow(sqrt(pow(x_gl(ii, jj) - xC(i, j), 2) + pow(y_gl(ii, jj) - yC(i, j), 2)), 3)
                                                      * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                      * (xi_1(i-1) + 1)/2 * (xi_2(j+1) - xi_2(j-1))/2;
                    }
                    if (i < nx-1)
                    {
                        std::vector<fastgl::QuadPair> gl_x(nx-i+2), gl_y(4);
                        arma::vec x1_gl(nx-i+2), x2_gl(4);
                        for (size_t ii = 0; ii < nx-i+2; ii++)
                        {
                            gl_x[ii] = fastgl::GLPair(nx-i+2, ii+1);
                            x1_gl(ii) = (1 - gl_x[ii].x()) * (1 - xi_1(i+1))/2 + xi_1(i+1);
                        }
                        for (size_t jj = 0; jj < 4; jj++)
                        {
                            gl_y[jj] = fastgl::GLPair(4, jj+1);
                            x2_gl(jj) = (1 - gl_y[jj].x()) * (xi_2(j+1) - xi_2(j-1))/2 + xi_2(j-1);
                        }
                        auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                        auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                        arma::mat J_gl = dx_gldx1%dy_gldx2 - dx_gldx2%dy_gldx1;
                        for (size_t ii = 0; ii < nx-i+2; ii++)
                            for (size_t jj = 0; jj < 4; jj++)
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        A(k1, p+q*nx) += gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                      / pow(sqrt(pow(x_gl(ii, jj) - xC(i, j), 2) + pow(y_gl(ii, jj) - yC(i, j), 2)), 3)
                                                      * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                      * (1 - xi_1(i+1))/2 * (xi_2(j+1) - xi_2(j-1))/2;
                    }
                }
            break;
        }
        case Analysis::nonlinear:
        {
            size_t m_rho = 2;
            break;
        }
    }
}