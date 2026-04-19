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

void Wing::regularIntegralLinear(size_t k, double xC, double yC, size_t nx_ii, size_t ny_jj, double xm, double xp, double ym, double yp)
{
    std::vector<fastgl::QuadPair> gl_x(nx_ii), gl_y(ny_jj);
    arma::vec x1_gl(nx_ii), x2_gl(ny_jj);
    for (size_t ii = 0; ii < nx_ii; ii++)
    {
        gl_x[ii] = fastgl::GLPair(nx_ii, ii+1);
        x1_gl(ii) = (1 - gl_x[ii].x()) * (xp - xm)/2 + xm;
    }
    for (size_t jj = 0; jj < ny_jj; jj++)
    {
        gl_y[jj] = fastgl::GLPair(ny_jj, jj+1);
        x2_gl(jj) = (1 - gl_y[jj].x()) * (yp - ym)/2 + ym;
    }
    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
    for (size_t ii = 0; ii < nx_ii; ii++)
        for (size_t jj = 0; jj < ny_jj; jj++)
        {
            arma::vec::fixed<2> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC};
            double r3 = pow(norm(r), 3);
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double T2     = boost::math::chebyshev_t(q, x2_gl(jj));
                double dT2dx2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double T1     = boost::math::chebyshev_t(p, x1_gl(ii));
                    double dT1dx1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                    double dmudx1 = dT1dx1 * T2;
                    double dmudx2 = T1 * dT2dx2;
                    A(k, p+q*nx) += gl_x[ii].weight*gl_y[jj].weight/r3*(xp - xm)*(yp - ym)/4
                                 *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                 - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                }
            }
        }
}

void Wing::regularIntegralNonlinear(size_t k, double xC, double yC, double zC, size_t nx_ii, size_t ny_jj, double xm, double xp, double ym, double yp)
{
    std::vector<fastgl::QuadPair> gl_x(nx_ii), gl_y(ny_jj);
    arma::vec x1_gl(nx_ii), x2_gl(ny_jj);
    for (size_t ii = 0; ii < nx_ii; ii++)
    {
        gl_x[ii] = fastgl::GLPair(nx_ii, ii+1);
        x1_gl(ii) = (1 - gl_x[ii].x()) * (xp - xm)/2 + xm;
    }
    for (size_t jj = 0; jj < ny_jj; jj++)
    {
        gl_y[jj] = fastgl::GLPair(ny_jj, jj+1);
        x2_gl(jj) = (1 - gl_y[jj].x()) * (yp - ym)/2 + ym;
    }
    arma::mat Tx = Lagrange::interpolationMatrix(x1, x1_gl);
    arma::mat Ty = Lagrange::interpolationMatrix(x2, x2_gl);
    arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, z, x1_gl, x2_gl);
    arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
    arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
    arma::mat dzdx1_gl = D1_gl*z_gl;
    arma::mat dzdx2_gl = z_gl*D2_gl.t();
    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
    auto [dxdx1_gl, dxdx2_gl, dydx1_gl, dydx2_gl] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
    arma::field<arma::mat> J_gl = {{dxdx1_gl, dxdx2_gl}, {dydx1_gl, dydx2_gl}};
    arma::cube e_c_gl = MetricCo(J_gl);
    arma::cube ec_gl  = MetricContra(e_c_gl);
    arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
    arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
    for (size_t ii = 0; ii < nx_ii; ii++)
        for (size_t jj = 0; jj < ny_jj; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC, z_gl(ii, jj) - zC};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dydx1_gl(ii, jj)*dzdx2_gl(ii, jj)-dzdx1_gl(ii, jj)*dydx2_gl(ii, jj),
                                                            dzdx1_gl(ii, jj)*dxdx2_gl(ii, jj)-dxdx1_gl(ii, jj)*dzdx2_gl(ii, jj),
                                                            dxdx1_gl(ii, jj)*dydx2_gl(ii, jj)-dydx1_gl(ii, jj)*dxdx2_gl(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dydx2_gl(ii, jj)*n_gl(2) - n_gl(1)*dzdx2_gl(ii, jj),-(dydx1_gl(ii, jj)*n_gl(2) - n_gl(1)*dzdx1_gl(ii, jj))},
                                            {-(dxdx2_gl(ii, jj)*n_gl(2) - n_gl(0)*dzdx2_gl(ii, jj)),dxdx1_gl(ii, jj)*n_gl(2) - n_gl(0)*dzdx1_gl(ii, jj)},
                                            {dxdx2_gl(ii, jj)*n_gl(1) - n_gl(0)*dydx2_gl(ii, jj),-(dxdx1_gl(ii, jj)*n_gl(1) - n_gl(0)*dydx1_gl(ii, jj))}};
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double T_1 = boost::math::chebyshev_t_prime(p, x1_gl(ii)) * boost::math::chebyshev_t(q, x2_gl(jj));
                    double T_2 = boost::math::chebyshev_t(p, x1_gl(ii)) * boost::math::chebyshev_t_prime(q, x2_gl(jj));
                    arma::vec::fixed<2> dmudxi = {T_1, T_2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3 * (xp - xm)/2 * (yp - ym)/2;
                    A(k, p+q*nx) += dot(q_mu, nC.row(k));
                }
        }
}

/**
 * @brief 
 * 
 */
void Wing::aerodynamicMatrix()
{
    arma::vec theta(n_theta);
    std::vector<fastgl::QuadPair> gl_theta(n_theta);
    for (size_t n = 0; n < n_theta; n++)
    {
        gl_theta[n] = fastgl::GLPair(n_theta, n+1);
        theta(n) =-arma::datum::pi * gl_theta[n].x();
    }
    std::vector<fastgl::QuadPair> gl_rho(m_rho);
    for (size_t m = 0; m < m_rho; m++)
        gl_rho[m] = fastgl::GLPair(m_rho, m+1);
    arma::vec cT = cos(theta);
    arma::vec sT = sin(theta);
    arma::vec scT = sin(theta)%cos(theta);
    arma::vec c2T = pow(cos(theta), 2)/2;
    arma::vec s2T = pow(sin(theta), 2)/2;
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
                    
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double T2 = boost::math::chebyshev_t(q, eta_2);
                        double dT2dx2  = boost::math::chebyshev_t_prime(q, eta_2);
                        double d2T2dx2dx2 = q<2 ? 0 : q/4.*((q+1)*boost::math::chebyshev_t(q-2, eta_2) - 2*q*T2 + (q-1)*boost::math::chebyshev_t(q+2, eta_2))/pow(1-pow(eta_2, 2), 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1 = boost::math::chebyshev_t(p, eta_1);
                            double dT1dx1  = boost::math::chebyshev_t_prime(p, eta_1);
                            double d2T1dx1dx1 = p<2 ? 0 : p/4.*((p+1)*boost::math::chebyshev_t(p-2, eta_1) - 2*p*T1 + (p-1)*boost::math::chebyshev_t(p+2, eta_1))/pow(1-pow(eta_1, 2), 2);
                            double T_1 = dT1dx1 * T2;
                            double T_2 = T1 * dT2dx2;
                            double T_11 = d2T1dx1dx1 * T2;
                            double T_12 = dT1dx1 * dT2dx2;
                            double T_22 = T1 * d2T2dx2dx2;
                            double gammax =-(dxdxi_2(i, j)*T_1 - dxdxi_1(i, j)*T_2);
                            double gammay =-(dydxi_2(i, j)*T_1 - dydxi_1(i, j)*T_2);
                            double dgammaxdxi_1 =-(d2xdxi_1dxi_2(i, j)*T_1 - d2xdxi_12(i, j)*T_2 + dxdxi_2(i, j)*T_11 - dxdxi_1(i, j)*T_12);
                            double dgammaxdxi_2 =-(d2xdxi_22(i, j)*T_1 - d2xdxi_1dxi_2(i, j)*T_2 + dxdxi_2(i, j)*T_12 - dxdxi_1(i, j)*T_22);
                            double dgammaydxi_1 =-(d2ydxi_1dxi_2(i, j)*T_1 - d2ydxi_12(i, j)*T_2 + dydxi_2(i, j)*T_11 - dydxi_1(i, j)*T_12);
                            double dgammaydxi_2 =-(d2ydxi_22(i, j)*T_1 - d2ydxi_1dxi_2(i, j)*T_2 + dydxi_2(i, j)*T_12 - dydxi_1(i, j)*T_22);
                            arma::vec gamma_1x = dgammaxdxi_1*cT + dgammaxdxi_2*sT;
                            arma::vec gamma_1y = dgammaydxi_1*cT + dgammaydxi_2*sT;
                            double I0  = 0;
                            double I_1 = 0;
                            for (size_t n = 0; n < n_theta; n++)
                            {
                                arma::vec::fixed<2> A_t = dXdxi_1*cT(n) + dXdxi_2*sT(n);
                                arma::vec::fixed<2> B_t = d2Xdxi_12*c2T(n) + d2Xdxi_1dxi_2*scT(n) + d2Xdxi_22*s2T(n);
                                double C_t = dot(A_t, B_t);
                                double A_abs = norm(A_t);
                                double A3 = pow(A_abs, 3);
                                double A5 = pow(A_abs, 5);
                                arma::vec::fixed<2> c1 = A_t/A3;
                                arma::vec::fixed<2> c2 = B_t/A3 - 3*A_t*C_t/A5;
                                arma::vec::fixed<2> c3 =-3*B_t*C_t/A5;
                                double F_1 = gammax*c1(1) - gammay*c1(0);
                                double F0  = gammax*c2(1) - gammay*c2(0) + gamma_1x(n)*c1(1) - gamma_1y(n)*c1(0);
                                double F1  = gammax*c3(1) - gammay*c3(0) + gamma_1x(n)*c2(1) - gamma_1y(n)*c2(0);
                                double F2  =                               gamma_1x(n)*c3(1) - gamma_1y(n)*c3(0);

                                double I_r = 0;
                                for (size_t m = 0; m < m_rho; m++)
                                {
                                    double rho = rho_tilde(n)/2*(1 - gl_rho[m].x());
                                    I_r += gl_rho[m].weight * (F0 + F1*rho + F2*pow(rho, 2));
                                }
                                I_r *= rho_tilde(n)/2;
                                I0  += gl_theta[n].weight * I_r;
                                I_1 += gl_theta[n].weight * F_1*log(fabs(rho_tilde(n)*A_abs));
                            }
                            A(k1, p+q*nx) = arma::datum::pi*(I0 + I_1);
                        }
                    }
                }
            }
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;
                    size_t k1 = i + j*nx;
                    if (j > 1)
                        regularIntegralLinear(k1, xC(i, j), yC(i, j), i+2, j+2, -1, x_right, -1, y_lower);
                    if (j < ny-1)
                        regularIntegralLinear(k1, xC(i, j), yC(i, j), nx+2-i, ny+2-j, x_left, 1, y_upper, 1);
                    if (i > 1)
                        regularIntegralLinear(k1, xC(i, j), yC(i, j), i+2, ny+2-j, -1, x_left, y_lower, 1);
                    if (i < nx-1)
                        regularIntegralLinear(k1, xC(i, j), yC(i, j), nx+2-i, j+2, x_right, 1, -1, y_upper);
                }
            }
            break;
        }
        case Analysis::nonlinear:
        {
            arma::mat Tx = Lagrange::interpolationMatrix(x1, xi_1);
            arma::mat Ty = Lagrange::interpolationMatrix(x2, xi_2);
            arma::mat zC = Lagrange::interpolation2D(Tx, Ty, z, xi_1, xi_2);
            arma::mat dzdxi_1 = D1*zC;
            arma::mat dzdxi_2 = zC*D2.t();
            arma::mat d2zdxi_12 = D1*dzdxi_1;
            arma::mat d2zdxi_1dxi_2 = D1*dzdxi_2;
            arma::mat d2zdxi_22 = dzdxi_2*D2.t();

            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    double e_11 = e_c(i, j, 0);
                    double e_12 = e_c(i, j, 1);
                    double e_22 = e_c(i, j, 2);
                    double e11  =  ec(i, j, 0);
                    double e12  =  ec(i, j, 1);
                    double e22  =  ec(i, j, 2);
                    double e = e_11*e_22 - pow(e_12, 2);
                    double sqrt_a = sqrt(e*(1 + e11*pow(dzdxi_1(i, j), 2) + 2*e12*dzdxi_1(i, j)*dzdxi_2(i, j) + e22*pow(dzdxi_2(i, j), 2)));
                    nC.row(i+j*nx) = arma::rowvec::fixed<3>({dydxi_1(i, j)*dzdxi_2(i, j)-dzdxi_1(i, j)*dydxi_2(i, j),
                                                             dzdxi_1(i, j)*dxdxi_2(i, j)-dxdxi_1(i, j)*dzdxi_2(i, j),
                                                             dxdxi_1(i, j)*dydxi_2(i, j)-dydxi_1(i, j)*dxdxi_2(i, j)})/sqrt_a;
                }
            arma::mat dnxdxi_1 = D1*reshape(nC.col(0), nx, ny);
            arma::mat dnxdxi_2 = reshape(nC.col(0), nx, ny)*D2.t();
            arma::mat dnydxi_1 = D1*reshape(nC.col(1), nx, ny);
            arma::mat dnydxi_2 = reshape(nC.col(1), nx, ny)*D2.t();
            arma::mat dnzdxi_1 = D1*reshape(nC.col(2), nx, ny);
            arma::mat dnzdxi_2 = reshape(nC.col(2), nx, ny)*D2.t();
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double eta_2 = xi_2(j);

                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;

                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    size_t k1 = i+j*nx;
                    double eta_1 = xi_1(i);
                    arma::vec::fixed<3> dXdxi_1 = {dxdxi_1(i, j), dydxi_1(i, j), dzdxi_1(i, j)};
                    arma::vec::fixed<3> dXdxi_2 = {dxdxi_2(i, j), dydxi_2(i, j), dzdxi_2(i, j)};
                    arma::vec::fixed<3> d2Xdxi_12 = {d2xdxi_12(i, j), d2ydxi_12(i, j), d2zdxi_12(i, j)};
                    arma::vec::fixed<3> d2Xdxi_1dxi_2 = {d2xdxi_1dxi_2(i, j), d2ydxi_1dxi_2(i, j), d2zdxi_1dxi_2(i, j)};
                    arma::vec::fixed<3> d2Xdxi_22 = {d2xdxi_22(i, j), d2ydxi_22(i, j), d2zdxi_22(i, j)};

                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, eta_1, eta_2, theta); 

                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double T2 = boost::math::chebyshev_t(q, eta_2);
                        double dT2dx2  = boost::math::chebyshev_t_prime(q, eta_2);
                        double d2T2dx2dx2 = q<2? 0 : q/4.*((q+1)*boost::math::chebyshev_t(q-2, eta_2) - 2*q*T2 + (q-1)*boost::math::chebyshev_t(q+2, eta_2))/pow(1-pow(eta_2, 2), 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1 = boost::math::chebyshev_t(p, eta_1);
                            double dT1dx1  = boost::math::chebyshev_t_prime(p, eta_1);
                            double d2T1dx1dx1 = p<2 ? 0 : p/4.*((p+1)*boost::math::chebyshev_t(p-2, eta_1) - 2*p*T1 + (p-1)*boost::math::chebyshev_t(p+2, eta_1))/pow(1-pow(eta_1, 2), 2);
                            double T_1 = dT1dx1 * T2;
                            double T_2 = T1 * dT2dx2;
                            double T_11 = d2T1dx1dx1 * T2;
                            double T_12 = dT1dx1 * dT2dx2;
                            double T_22 = T1 * d2T2dx2dx2;
                            double dmudx = (dydxi_2(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_2(i, j))*T_1
                                         - (dydxi_1(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_1(i, j))*T_2;
                            double dmudy = (dxdxi_1(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_1(i, j))*T_2
                                         - (dxdxi_2(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_2(i, j))*T_1;
                            double dmudz = (dxdxi_2(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_2(i, j))*T_1
                                         - (dxdxi_1(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_1(i, j))*T_2;
                            arma::vec::fixed<3> gamma_0 = {dmudy*nC(k1, 2) - dmudz*nC(k1, 1),
                                                           dmudz*nC(k1, 0) - dmudx*nC(k1, 2),
                                                           dmudx*nC(k1, 1) - dmudy*nC(k1, 0)};
                            double dmudxdxi_1 = (dydxi_2(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_2(i, j))*T_11
                                              + (d2ydxi_1dxi_2(i, j)*nC(k1, 2) + dydxi_2(i, j)*dnzdxi_1(i, j) - dnydxi_1(i, j)*dzdxi_2(i, j) - nC(k1, 1)*d2zdxi_1dxi_2(i, j))*T_1
                                              - (dydxi_1(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_1(i, j))*T_12
                                              - (d2ydxi_12(i, j)*nC(k1, 2) + dydxi_1(i, j)*dnzdxi_1(i, j) - dnydxi_1(i, j)*dzdxi_1(i, j) - nC(k1, 1)*d2zdxi_12(i, j))*T_2;
                            double dmudxdxi_2 = (dydxi_2(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_2(i, j))*T_12
                                              + (d2ydxi_22(i, j)*nC(k1, 2) + dydxi_2(i, j)*dnzdxi_2(i, j) - dnydxi_2(i, j)*dzdxi_2(i, j) - nC(k1, 1)*d2zdxi_22(i, j))*T_1
                                              - (dydxi_1(i, j)*nC(k1, 2) - nC(k1, 1)*dzdxi_1(i, j))*T_22
                                              - (d2ydxi_1dxi_2(i, j)*nC(k1, 2) + dydxi_1(i, j)*dnzdxi_2(i, j) - dnydxi_2(i, j)*dzdxi_1(i, j) - nC(k1, 1)*d2zdxi_1dxi_2(i, j))*T_2;
                            double dmudydxi_1 = (dxdxi_1(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_1(i, j))*T_12
                                              + (d2xdxi_12(i, j)*nC(k1, 2) + dxdxi_1(i, j)*dnzdxi_1(i, j) - dnxdxi_1(i, j)*dzdxi_1(i, j) - nC(k1, 0)*d2zdxi_12(i, j))*T_2
                                              - (dxdxi_2(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_2(i, j))*T_11
                                              - (d2xdxi_1dxi_2(i, j)*nC(k1, 2) + dxdxi_2(i, j)*dnzdxi_1(i, j) - dnxdxi_1(i, j)*dzdxi_2(i, j) - nC(k1, 0)*d2zdxi_1dxi_2(i, j))*T_1;
                            double dmudydxi_2 = (dxdxi_1(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_1(i, j))*T_22
                                              + (d2xdxi_1dxi_2(i, j)*nC(k1, 2) + dxdxi_1(i, j)*dnzdxi_2(i, j) - dnxdxi_2(i, j)*dzdxi_1(i, j) - nC(k1, 0)*d2zdxi_1dxi_2(i, j))*T_2
                                              - (dxdxi_2(i, j)*nC(k1, 2) - nC(k1, 0)*dzdxi_2(i, j))*T_12
                                              - (d2xdxi_22(i, j)*nC(k1, 2) + dxdxi_2(i, j)*dnzdxi_2(i, j) - dnxdxi_2(i, j)*dzdxi_2(i, j) - nC(k1, 0)*d2zdxi_22(i, j))*T_1;
                            double dmudzdxi_1 = (dxdxi_2(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_2(i, j))*T_11
                                              + (d2xdxi_1dxi_2(i, j)*nC(k1, 1) + dxdxi_2(i, j)*dnydxi_1(i, j) - dnxdxi_1(i, j)*dydxi_2(i, j) - nC(k1, 0)*d2ydxi_1dxi_2(i, j))*T_1
                                              - (dxdxi_1(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_1(i, j))*T_12
                                              - (d2xdxi_12(i, j)*nC(k1, 1) + dxdxi_1(i, j)*dnydxi_1(i, j) - dnxdxi_1(i, j)*dydxi_1(i, j) - nC(k1, 0)*d2ydxi_12(i, j))*T_2;
                            double dmudzdxi_2 = (dxdxi_2(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_2(i, j))*T_12
                                              + (d2xdxi_22(i, j)*nC(k1, 1) + dxdxi_2(i, j)*dnydxi_2(i, j) - dnxdxi_2(i, j)*dydxi_2(i, j) - nC(k1, 0)*d2ydxi_22(i, j))*T_1
                                              - (dxdxi_1(i, j)*nC(k1, 1) - nC(k1, 0)*dydxi_1(i, j))*T_22
                                              - (d2xdxi_1dxi_2(i, j)*nC(k1, 1) + dxdxi_1(i, j)*dnydxi_2(i, j) - dnxdxi_2(i, j)*dydxi_1(i, j) - nC(k1, 0)*d2ydxi_1dxi_2(i, j))*T_2;
                            arma::vec::fixed<3> gamma_1_1 = {dmudydxi_1*nC(k1, 2) + dmudy*dnzdxi_1(i, j) - dmudzdxi_1*nC(k1, 1) - dmudz*dnydxi_1(i, j),
                                                             dmudzdxi_1*nC(k1, 0) + dmudz*dnxdxi_1(i, j) - dmudxdxi_1*nC(k1, 2) - dmudx*dnzdxi_1(i, j),
                                                             dmudxdxi_1*nC(k1, 1) + dmudx*dnydxi_1(i, j) - dmudydxi_1*nC(k1, 0) - dmudy*dnxdxi_1(i, j)};
                            arma::vec::fixed<3> gamma_1_2 = {dmudydxi_2*nC(k1, 2) + dmudy*dnzdxi_2(i, j) - dmudzdxi_2*nC(k1, 1) - dmudz*dnydxi_2(i, j),
                                                             dmudzdxi_2*nC(k1, 0) + dmudz*dnxdxi_2(i, j) - dmudxdxi_2*nC(k1, 2) - dmudx*dnzdxi_2(i, j),
                                                             dmudxdxi_2*nC(k1, 1) + dmudx*dnydxi_2(i, j) - dmudydxi_2*nC(k1, 0) - dmudy*dnxdxi_2(i, j)};
                            arma::mat gamma_1 = gamma_1_1*cT.t() + gamma_1_2*sT.t();
                            arma::vec::fixed<3> I0(arma::fill::zeros);
                            arma::vec::fixed<3> I_1(arma::fill::zeros);
                            arma::vec mu1 = T_1*cT + T_2*sT;
                            for (size_t n = 0; n < n_theta; n++)
                            {
                                arma::vec::fixed<3> A_t = dXdxi_1*cT(n) + dXdxi_2*sT(n);
                                arma::vec::fixed<3> B_t = d2Xdxi_12*c2T(n) + d2Xdxi_1dxi_2*scT(n) + d2Xdxi_22*s2T(n);
                                double C_t  = dot(A_t, B_t);
                                double A_abs = norm(A_t);
                                double A3 = pow(A_abs, 3);
                                double A5 = pow(A_abs, 5);
                                arma::vec::fixed<3> c1 = A_t/A3;
                                arma::vec::fixed<3> c2 = B_t/A3 - 3*A_t*C_t/A5;
                                arma::vec::fixed<3> c3 = -3*B_t*C_t/A5;
                                arma::vec::fixed<3> F_1 = cross(gamma_0, c1);
                                arma::vec::fixed<3> F0  = cross(gamma_0, c2) + cross(gamma_1.col(n), c1);
                                arma::vec::fixed<3> F1  = cross(gamma_0, c3) + cross(gamma_1.col(n), c2);
                                arma::vec::fixed<3> F2  =                      cross(gamma_1.col(n), c3);

                                arma::vec::fixed<3> I_r(arma::fill::zeros);
                                for (size_t m = 0; m < m_rho; m++)
                                {
                                    double rho = rho_tilde(n)*(1 - gl_rho[m].x())/2;
                                    I_r += gl_rho[m].weight * (F0 + F1*rho + F2*pow(rho, 2));
                                }
                                I_r *= rho_tilde(n)/2;
                                I0  += gl_theta[n].weight * I_r;
                                I_1 += gl_theta[n].weight * F_1*log(fabs(rho_tilde(n)*A_abs));
                            }
                            arma::vec::fixed<3> q_mu = arma::datum::pi*(I0 + I_1);
                            A(k1, p+q*nx) = dot(q_mu, nC.row(k1));
                        }
                    }
                }
            }
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double y_lower = j>1    ? xi_2(j-1) :-1;
                double y_upper = j<ny-2 ? xi_2(j+1) : 1;
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double x_left  = i>1    ? xi_1(i-1) :-1;
                    double x_right = i<nx-2 ? xi_1(i+1) : 1;
                    size_t k1 = i + j*nx;
                    if (j > 1)
                        regularIntegralNonlinear(k1, xC(i, j), yC(i, j), zC(i, j), i+2, j+2, -1, x_right, -1, y_lower);
                    if (j < ny-1)
                        regularIntegralNonlinear(k1, xC(i, j), yC(i, j), zC(i, j), nx+2-i, ny+2-j, x_left, 1, y_upper, 1);
                    if (i > 1)
                        regularIntegralNonlinear(k1, xC(i, j), yC(i, j), zC(i, j), i+2, ny+2-j, -1, x_left, y_lower, 1);
                    if (i < nx-1)
                        regularIntegralNonlinear(k1, xC(i, j), yC(i, j), zC(i, j), nx+2-i, j+2, x_right, 1, -1, y_upper);
                }
            }
            break;
        }
    }
}