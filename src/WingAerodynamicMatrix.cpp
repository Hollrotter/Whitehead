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
    auto [xC, yC] = Lagrange::TransfiniteQuadMap(xi_1, xi_2, chi);
    auto [dxdxi_1, dxdxi_2, dydxi_1, dydxi_2] = Lagrange::TransfiniteQuadMetrics(xi_1, xi_2, chi);
    arma::mat d2xdxi_12     = D1*dxdxi_1;
    arma::mat d2xdxi_1dxi_2 = dxdxi_1*D2.t();
    arma::mat d2xdxi_22     = dxdxi_2*D2.t();
    arma::mat d2ydxi_12     = D1*dydxi_1;
    arma::mat d2ydxi_1dxi_2 = D1*dydxi_2;
    arma::mat d2ydxi_22     = dydxi_2*D2.t();
    switch (analysis)
    {
        case Analysis::linear:
        {
            size_t m_rho = 1;
            std::vector<fastgl::QuadPair> gl_rho(m_rho);
            for (size_t m = 0; m < m_rho; m++)
                gl_rho[m] = fastgl::GLPair(m_rho, m+1);
            arma::mat detJ = dxdxi_1%dydxi_2 - dxdxi_2%dydxi_1;
            arma::mat dJdxi_1 = D1*detJ;
            arma::mat dJdxi_2 = detJ*D2.t();
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double eta_2 = xi_2(j);

                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double eta_1 = xi_1(i);

                    arma::vec::fixed<2> dXdxi_1 = {dxdxi_1(i, j), dydxi_1(i, j)};
                    arma::vec::fixed<2> dXdxi_2 = {dxdxi_2(i, j), dydxi_2(i, j)};

                    arma::vec::fixed<2> d2Xdxi_12     = {d2xdxi_12(i, j),     d2ydxi_12(i, j)};
                    arma::vec::fixed<2> d2Xdxi_1dxi_2 = {d2xdxi_1dxi_2(i, j), d2ydxi_1dxi_2(i, j)};
                    arma::vec::fixed<2> d2Xdxi_22     = {d2xdxi_22(i, j),     d2ydxi_22(i, j)};

                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, eta_1, eta_2, theta);

                    size_t k1 = i+j*nx;
                    b.row(k1) =-4*arma::datum::pi*alpha; // Must be corrected later!

                    double J0 = detJ(i, j);
                    arma::vec J1 = dJdxi_1(i, j)*cT + dJdxi_2(i, j)*sT;
                    
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T_1 = boost::math::chebyshev_t_prime(p, eta_1) * boost::math::chebyshev_t(q, eta_2);
                            double T_2 = boost::math::chebyshev_t(p, eta_1)       * boost::math::chebyshev_t_prime(q, eta_2);
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
            std::vector<fastgl::QuadPair> gl_rho(m_rho);
            for (size_t m = 0; m < m_rho; m++)
                gl_rho[m] = fastgl::GLPair(m_rho, m+1);
            arma::mat detJ;
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double eta_2 = xi_2(j);

                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;

                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double eta_1 = xi_1(i);
                    arma::vec::fixed<3> dXdxi_1;
                    arma::vec::fixed<3> dXdxi_2;
                    arma::vec::fixed<3> d2Xdxi_12;
                    arma::vec::fixed<3> d2Xdxi_1dxi_2;
                    arma::vec::fixed<3> d2Xdxi_22;

                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, eta_1, eta_2, theta);

                    size_t k1 = i+j*nx;

                    arma::vec J0 = detJ(i, j)*nC;
                    double dJdxi_1;
                    double dJdxi_2;
                    arma::vec J1 = dJdxi_1*cT + dJdxi_2*sT;
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T_1 = boost::math::chebyshev_t_prime(p, eta_1) * boost::math::chebyshev_t(q, eta_2);
                            double T_2 = boost::math::chebyshev_t(p, eta_1)       * boost::math::chebyshev_t_prime(q, eta_2);
                            double mu0 = boost::math::chebyshev_t(p, eta_1) * boost::math::chebyshev_t(q, eta_2);
                            arma::vec::fixed<3> I0(arma::fill::zeros);
                            arma::vec::fixed<3> I_1(arma::fill::zeros);
                            arma::vec::fixed<3> I_2(arma::fill::zeros);
                            arma::vec mu1 = T_1*cT + T_2*sT;
                            for (size_t n = 0; n < n_theta; n++)
                            {
                                arma::vec::fixed<3> A_t = dXdxi_1*cT(n) + dXdxi_2*sT(n);
                                arma::vec::fixed<3> B_t = d2Xdxi_12*c2T(n) + d2Xdxi_1dxi_2*scT(n) + d2Xdxi_22*s2T(n);
                                double C_t  = dot(A_t, B_t);
                                double A_J1 = dot(A_t, J1);
                                double B_J0 = dot(B_t, J0);
                                double B_J1 = dot(B_t, J1);
                                double A_abs = norm(A_t);
                                double A2 = pow(A_abs, 2);
                                double A3 = pow(A_abs, 3);
                                double A4 = pow(A2, 2);
                                double A5 = A2*A3;
                                double A7 = A2*A5;
                                arma::vec::fixed<3> F_2 = J0*mu0/A3;
                                arma::vec::fixed<3> F_1 = ((J1 - 3/A2*(C_t*J0 + (A_J1 + B_J0)*A_t))*mu0 + J0*mu1(n))/A3;
                                arma::vec::fixed<3> F0 =-3/A5*(((A_J1 + B_J0)*B_t + A_t*B_J1 + C_t*J1 - A_t*(A_J1 + B_J0)*5*C_t/A2)*mu0
                                                       + (C_t*J0 + A_t*(A_J1 + B_J0) - J1/3*A2)*mu1(n));
                                arma::vec::fixed<3> F1 =-3/A5*((B_t*B_J1 - 5*C_t/A2*((A_J1 + B_J0)*B_t + A_t*B_J1))*mu0
                                                       + (J1*C_t - 5*A_t*(A_J1 + B_J0)*C_t/A2 + (A_J1 + B_J0)*B_t + A_t*B_J1)*mu1(n));
                                arma::vec::fixed<3> F2 = 3/A5*(B_J1*5*C_t/A2*B_t*mu0
                                                       + (((A_J1 + B_J0)*B_t + A_t*B_J1)*5*C_t/A2- B_J1*B_t)*mu1(n));
                                arma::vec::fixed<3> F3 = 15*B_t*B_J1*C_t/A7*mu1(n);

                                arma::vec::fixed<3> I_r(arma::fill::zeros);
                                for (size_t m = 0; m < m_rho; m++)
                                {
                                    double rho = rho_tilde(n)*(1 - gl_rho[m].x())/2;
                                    I_r += gl_rho[m].weight * (F0 + F1*rho + F2*pow(rho, 2) + F3*pow(rho, 3));
                                }
                                I_r *= rho_tilde(n)/2;
                                I0  += gl_theta[n].weight * I_r;
                                I_1 += gl_theta[n].weight * F_1*log(fabs(rho_tilde(n)*A_abs));
                                I_2 -= gl_theta[n].weight * F_2*(1/rho_tilde(n) + C_t/A2);
                            }
                            arma::vec::fixed<3> q_mu = arma::datum::pi*(I0 + I_1 + I_2);
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
                        arma::mat J_gl;
                        double zC, z_gl;
                        for (size_t ii = 0; ii < nx; ii++)
                            for (size_t jj = 0; jj < j+2; jj++)
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl - zC};
                                double r_abs = norm(r);
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                       * (nC - 3*(dot(r, nC))*r/pow(r_abs, 2))/pow(r_abs, 3)
                                                       * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                       * (xi_2(j-1) + 1)/2;
                            }
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
                        arma::mat J_gl;
                        double zC, z_gl;
                        for (size_t ii = 0; ii < nx; ii++)
                            for (size_t jj = 0; jj < ny-j+2; jj++)
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl - zC};
                                double r_abs = norm(r);
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                       * (nC - 3*(dot(r, nC))*r/pow(r_abs, 2))/pow(r_abs, 3)
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
                        arma::mat J_gl;
                        double zC, z_gl;
                        for (size_t ii = 0; ii < i+2; ii++)
                            for (size_t jj = 0; jj < 4; jj++)
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl - zC};
                                double r_abs = norm(r);
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                       * (nC - 3*(dot(r, nC))*r/pow(r_abs, 2))/pow(r_abs, 3)
                                                       * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                       * (xi_1(i-1) + 1)/2 * (xi_2(j+1) - xi_2(j-1))/2;
                            }
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
                        arma::mat J_gl;
                        double zC, z_gl;
                        for (size_t ii = 0; ii < nx-i+2; ii++)
                            for (size_t jj = 0; jj < 4; jj++)
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl - zC};
                                double r_abs = norm(r);
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * J_gl(ii, jj)
                                                       * (nC - 3*(dot(r, nC))*r/pow(r_abs, 2))/pow(r_abs, 3)
                                                       * boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj))
                                                       * (1 - xi_1(i+1))/2 * (xi_2(j+1) - xi_2(j-1))/2;
                            }
                    }
                }
            break;
        }
    }
}