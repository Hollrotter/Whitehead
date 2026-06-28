#include "Wing.hpp"

arma::vec Wing::externalContour(double x_min, double x_max, double y_min, double y_max, double eta_1, double eta_2, arma::vec theta)
{
    arma::vec rho(theta.size(), arma::fill::none);
    for (size_t n = 0; n < theta.size(); n++)
        if ((-arma::datum::pi < theta(n) && theta(n) <=-arma::datum::pi/2) || (arma::datum::pi < theta(n) && theta(n) <= 3*arma::datum::pi/2))
            rho(n) = std::min(fabs((x_min-eta_1)/cos(theta(n))), fabs((y_min-eta_2)/sin(theta(n))));
        else if ((-arma::datum::pi/2 < theta(n) && theta(n) <= 0) || (3*arma::datum::pi/2 < theta(n) && theta(n) <= 2*arma::datum::pi))
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
    arma::vec x1_gl(nx_ii, arma::fill::none), x2_gl(ny_jj, arma::fill::none);
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
    for (size_t jj = 0; jj < ny_jj; jj++)
    {
        double t2   = 1; // First Chebyshev Polynomial (j)
        double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

        double dt2   = 0;
        double dt2p1 = 1;
        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
        {
            for (size_t ii = 0; ii < nx_ii; ii++)
            {
                arma::vec::fixed<2> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC};
                double r3 = pow(norm(r), 3);
                double t1   = 1; // First Chebyshev Polynomial (i)
                double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = 0;
                double dt1p1 = 1;
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double dmudx1 = dt1 *  t2;
                    double dmudx2 =  t1 * dt2;
                    A(k, p+q*nx) += gl_x[ii].weight*gl_y[jj].weight/r3*(xp - xm)*(yp - ym)/4
                                    *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                    - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                    std::swap(t1, t1p1); // Swap order so the oldes Chebyshev Polynomial will be overwritten
                    t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1); // Calculate next Chebyshev Polynomial

                    std::swap(dt1, dt1p1);
                    dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                }
            }
            std::swap(t2, t2p1); // Swap order so the oldest Chebyshev Polynomial will be overwritten
            t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1); // Calculate next Chebyshev Polynomial

            std::swap(dt2, dt2p1);
            dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
        }
    }
}

void Wing::regularIntegralNonlinear(size_t k, double xC, double yC, double zC, size_t nx_ii, size_t ny_jj, double xm, double xp, double ym, double yp)
{
    std::vector<fastgl::QuadPair> gl_x(nx_ii), gl_y(ny_jj);
    arma::vec x1_gl(nx_ii, arma::fill::none), x2_gl(ny_jj, arma::fill::none);
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
    arma::mat Tx_gl = Lagrange::interpolationMatrix(x1, x1_gl);
    arma::mat Ty_gl = Lagrange::interpolationMatrix(x2, x2_gl);
    arma::mat z_gl  = Lagrange::interpolation2D(Tx_gl, Ty_gl, z, x1_gl, x2_gl);
    arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
    arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
    arma::mat dz_gldx1 = D1_gl*z_gl;
    arma::mat dz_gldx2 = z_gl*D2_gl.t();
    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
    arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
    arma::cube e_c_gl = MetricCo(J_gl);
    arma::cube ec_gl  = MetricContra(e_c_gl);
    arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
    arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
    for (size_t ii = 0; ii < nx_ii; ii++)
        for (size_t jj = 0; jj < ny_jj; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC, z_gl(ii, jj) - zC};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                            {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
            double t2   = 1; // First Chebyshev Polynomial (j)
            double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = 0;
            double dt2p1 = 1;
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double t1   = 1; // First Chebyshev Polynomial (i)
                double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = 0;
                double dt1p1 = 1;
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    arma::vec::fixed<2> dmudxi = {dt1 * t2, t1 * dt2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3 * (xp - xm) * (yp - ym)/4;
                    A(k, p+q*nx) += dot(q_mu, nC.row(k));
                    std::swap(t1, t1p1);
                    t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                }
                std::swap(t2, t2p1);
                t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
            }
        }
}

/**
 * @brief 
 * 
 */
void Wing::aerodynamicMatrix()
{
    arma::vec theta = arma::linspace(0, arma::datum::tau, n_theta);
    double dtheta = theta(1) - theta(0);
    arma::vec cT = cos(theta);
    arma::vec sT = sin(theta);
    arma::vec scT = sT%cT;
    arma::vec c2T = pow(cT, 2)/2;
    arma::vec s2T = pow(sT, 2)/2;
    auto [J11, J12, J21, J22] = Lagrange::TransfiniteQuadMetrics(xi_1, xi_2, chi);
    arma::mat dJ11dxi_1 = D1*J11;
    arma::mat dJ11dxi_2 = J11*D2.t();
    arma::mat dJ12dxi_1 = dJ11dxi_2;
    arma::mat dJ12dxi_2 = J12*D2.t();
    arma::mat dJ21dxi_1 = D1*J21;
    arma::mat dJ22dxi_1 = D1*J22;
    arma::mat dJ22dxi_2 = J22*D2.t();
    arma::mat dJ21dxi_2 = dJ22dxi_1;
    switch (analysis)
    {
        case Analysis::linear:
        {
            #pragma omp parallel for reduction(+:A)
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double y_lower = std::max(-1., xi_2(j)-delta/2);
                double y_upper = std::min(xi_2(j)+delta/2, 1.);
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double x_left  = std::max(-1., xi_1(i)-delta/2);
                    double x_right = std::min(xi_1(i)+delta/2, 1.);

                    arma::vec::fixed<2> dXdxi_1 = {J11(i, j), J21(i, j)};
                    arma::vec::fixed<2> dXdxi_2 = {J12(i, j), J22(i, j)};
                    arma::vec::fixed<2> d2Xdxi_12     = {dJ11dxi_1(i, j), dJ21dxi_1(i, j)};
                    arma::vec::fixed<2> d2Xdxi_1dxi_2 = {dJ11dxi_2(i, j), dJ22dxi_1(i, j)};
                    arma::vec::fixed<2> d2Xdxi_22     = {dJ12dxi_2(i, j), dJ22dxi_2(i, j)};

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, xi_1(i), xi_2(j), theta);

                    size_t k = i+j*nx;
                    double d2T2   = 0;
                    double d2T2p1 = 0;
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double d2T1   = 0;
                        double d2T1p1 = 0;
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double m1  = dT1(i, p) *  T2(j, q);
                            double m2  =  T1(i, p) * dT2(j, q);
                            double m11 =      d2T1 *  T2(j, q);
                            double m12 = dT1(i, p) * dT2(j, q);
                            double m22 =  T1(i, p) * d2T2;

                            double gammax = J11(i, j)*m2 - J12(i, j)*m1;
                            double gammay = J21(i, j)*m2 - J22(i, j)*m1;
                            double dgammaxdxi_1 = dJ11dxi_1(i, j)*m2 - dJ12dxi_1(i, j)*m1 - J12(i, j)*m11 + J11(i, j)*m12;
                            double dgammaxdxi_2 = dJ11dxi_2(i, j)*m2 - dJ12dxi_2(i, j)*m1 - J12(i, j)*m12 + J11(i, j)*m22;
                            double dgammaydxi_1 = dJ21dxi_1(i, j)*m2 - dJ22dxi_1(i, j)*m1 - J22(i, j)*m11 + J21(i, j)*m12;
                            double dgammaydxi_2 = dJ21dxi_2(i, j)*m2 - dJ22dxi_2(i, j)*m1 - J22(i, j)*m12 + J21(i, j)*m22;
                            arma::vec gamma_1x = dgammaxdxi_1*cT + dgammaxdxi_2*sT;
                            arma::vec gamma_1y = dgammaydxi_1*cT + dgammaydxi_2*sT;
                            double I0  = 0;
                            double I_1 = 0;
                            for (size_t n = 0; n < n_theta-1; n++)
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
                                I0  += F0*rho_tilde(n) + F1*pow(rho_tilde(n), 2)/2 + F2*pow(rho_tilde(n), 3)/3;
                                I_1 += F_1*log(fabs(rho_tilde(n)*A_abs));
                            }
                            A(k, p+q*nx) = dtheta*(I0 + I_1);
                            std::swap(d2T1, d2T1p1);
                            if (p == 0)
                                d2T1p1 = 4;
                            else if (p < nx-1)
                                d2T1p1 = (p+2)*(2*dT1(i, p+1) + d2T1p1/p);
                        }
                        std::swap(d2T2, d2T2p1);
                        if (q == 0)
                            d2T2p1 = 4;
                        else if (q < ny-1)
                            d2T2p1 = (q+2)*(2*dT2(j, q+1) + d2T2p1/q);
                    }
                    if (y_lower > -1)
                        regularIntegralLinear(k, xC(i, j), yC(i, j), std::ceil((x_right+1)/2*nx), std::ceil((y_lower+1)/2*ny), -1, x_right, -1, y_lower);
                    if (y_upper < 1)
                        regularIntegralLinear(k, xC(i, j), yC(i, j), std::ceil((1-x_left )/2*nx), std::ceil((1-y_upper)/2*ny),  x_left,  1, y_upper,  1);
                    if (x_left > -1)
                        regularIntegralLinear(k, xC(i, j), yC(i, j), std::ceil((x_left+1 )/2*nx), std::ceil((1-y_lower)/2*ny), -1,  x_left, y_lower,  1);
                    if (x_right < 1)
                        regularIntegralLinear(k, xC(i, j), yC(i, j), std::ceil((1-x_right)/2*nx), std::ceil((y_upper+1)/2*ny), x_right,  1, -1, y_upper);
                }
            }
            for (const Wake* w:wakes)
                for (size_t c = 0; c < 4; c++)
                    if (chi[c] == w->chi)
                    {
                        std::vector<fastgl::QuadPair> gl_x_w(nx), gl_y_w(ny);
                        arma::vec x1_gl_w(nx, arma::fill::none), x2_gl_w(ny, arma::fill::none);
                        for (size_t ii = 0; ii < nx; ii++)
                        {
                            gl_x_w[ii] = fastgl::GLPair(nx, ii+1);
                            x1_gl_w(ii) =-gl_x_w[ii].x();
                        }
                        for (size_t jj = 0; jj < ny; jj++)
                        {
                            gl_y_w[jj] = fastgl::GLPair(ny, jj+1);
                            x2_gl_w(jj) =-gl_y_w[jj].x();
                        }
                        arma::vec x1_gauss = Chebyshev::gaussLobatto(nx);
                        arma::vec x2_gauss = Chebyshev::gaussLobatto(ny);
                        arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                        arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                        if (c == 0) // South
                        {
                            arma::vec xw = T1_gauss * x.col(0);
                            arma::vec yw = T1_gauss * y.col(0);
                            for (size_t ii = 0; ii < nx; ii++)
                                for (size_t p = 0; p < nx; p++)
                                {
                                    double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                    for (size_t jj = 0; jj < ny; jj++)
                                    {
                                        double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                        for (size_t q = 0; q < ny; q++)
                                        {
                                            double t2 = pow(-1, q);
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    A(k, p+q*nx) += gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*t2*r(1)/pow(1 + x2_gl_w(jj), 2)/r3; 
                                                }
                                        }
                                    }
                                }
                        }
                        else if (c == 1) // East
                        {
                            arma::vec xw = T2_gauss * x.row(nx-1).t();
                            arma::vec yw = T2_gauss * y.row(nx-1).t();
                            for (size_t jj = 0; jj < ny; jj++)
                                for (size_t q = 0; q < ny; q++)
                                {
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                    for (size_t ii = 0; ii < nx; ii++)
                                    {
                                        double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                        for (size_t j = 1; j < ny-1; j++)
                                            for (size_t i = 1; i < nx-1; i++)
                                            {
                                                size_t k = i + j*nx;
                                                arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                double r3 = pow(norm(r), 3);
                                                double I = gl_x_w[ii].weight*gl_y_w[jj].weight*dt2*r(1)/pow(1 - x1_gl_w(ii), 2)/r3;
                                                for (size_t p = 0; p < nx; p++)
                                                    A(k, p+q*nx) += I;
                                            }
                                    }
                                }
                        }
                        else if (c == 2) // North
                        {
                            arma::vec xw = T1_gauss * x.col(ny-1);
                            arma::vec yw = T1_gauss * y.col(ny-1);
                            for (size_t ii = 0; ii < nx; ii++)
                                for (size_t p = 0; p < nx; p++)
                                {
                                    double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                    for (size_t jj = 0; jj < ny; jj++)
                                    {
                                        double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                        for (size_t j = 1; j < ny-1; j++)
                                            for (size_t i = 1; i < nx-1; i++)
                                            {
                                                size_t k = i + j*nx;
                                                arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                double r3 = pow(norm(r), 3);
                                                double I =-gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*r(1)/pow(1 - x2_gl_w(jj), 2)/r3;
                                                for (size_t q = 0; q < ny; q++)
                                                    A(k, p+q*nx) += I;
                                            }
                                    }
                                }
                        }
                        else // West
                        {
                            arma::vec xw = T2_gauss * x.row(0).t();
                            arma::vec yw = T2_gauss * y.row(0).t();
                            for (size_t jj = 0; jj < ny; jj++)
                                for (size_t q = 0; q < ny; q++)
                                {
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                    for (size_t ii = 0; ii < nx; ii++)
                                    {
                                        double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                        for (size_t p = 0; p < nx; p++)
                                        {
                                            double t1 = pow(-1, p);
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    A(k, p+q*nx) -= gl_x_w[ii].weight*gl_y_w[jj].weight*t1*dt2*r(1)/pow(1 + x1_gl_w(ii), 2)/r3;
                                                }
                                        }
                                    }
                                }
                        }
                    }
            if (sym == Symmetry::y)
            {
                std::vector<fastgl::QuadPair> gl_1(nx), gl_2(ny);
                arma::vec x1_gl(nx, arma::fill::none), x2_gl(ny, arma::fill::none);
                for (size_t ii = 0; ii < nx; ii++)
                {
                    gl_1[ii] = fastgl::GLPair(nx, ii+1);
                    x1_gl(ii) =-gl_1[ii].x();
                }
                for (size_t jj = 0; jj < ny; jj++)
                {
                    gl_2[jj] = fastgl::GLPair(ny, jj+1);
                    x2_gl(jj) =-gl_2[jj].x();
                }
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                y_gl     =-y_gl;
                dy_gldx1 =-dy_gldx1;
                dy_gldx2 =-dy_gldx2;
                #pragma omp parallel for reduction(+:A)
                for (size_t jj = 0; jj < ny; jj++)
                {
                    double t2   = 1; // First Chebyshev Polynomial (j)
                    double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

                    double dt2   = 0;
                    double dt2p1 = 1;
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        for (size_t ii = 0; ii < nx; ii++)
                        {
                            double t1   = 1; // First Chebyshev Polynomial (i)
                            double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                            double dt1   = 0;
                            double dt1p1 = 1;
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double dmudx1 = dt1 *  t2;
                                double dmudx2 =  t1 * dt2;
                                for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
                                    for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                                    {
                                        arma::vec::fixed<2> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j)};
                                        double r3 = pow(norm(r), 3);
                                        A(i+j*nx, p+q*nx) -= gl_1[ii].weight*gl_2[jj].weight/r3
                                                    *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                                    - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                                    }
                                std::swap(t1, t1p1); // Swap order so the oldes Chebyshev Polynomial will be overwritten
                                t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1); // Calculate next Chebyshev Polynomial

                                std::swap(dt1, dt1p1);
                                dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                            }
                        }
                        std::swap(t2, t2p1); // Swap order so the oldest Chebyshev Polynomial will be overwritten
                        t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1); // Calculate next Chebyshev Polynomial

                        std::swap(dt2, dt2p1);
                        dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
                    }
                }
                for (const Wake* w:wakes)
                    for (size_t c = 0; c < 4; c++)
                        if (chi[c] == w->chi)
                        {
                            std::vector<fastgl::QuadPair> gl_x(nx), gl_y(ny);
                            arma::vec x1_gl_w(nx, arma::fill::none), x2_gl_w(ny, arma::fill::none);
                            for (size_t ii = 0; ii < nx; ii++)
                            {
                                gl_x[ii] = fastgl::GLPair(nx, ii+1);
                                x1_gl_w(ii) =-gl_x[ii].x();
                            }
                            for (size_t jj = 0; jj < ny; jj++)
                            {
                                gl_y[jj] = fastgl::GLPair(ny, jj+1);
                                x2_gl_w(jj) =-gl_y[jj].x();
                            }
                            arma::vec x1_gauss = Chebyshev::gaussLobatto(nx);
                            arma::vec x2_gauss = Chebyshev::gaussLobatto(ny);
                            arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                            arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                            if (c == 0) // South
                            {
                                arma::vec xw = T1_gauss * x.col(0);
                                arma::vec yw =-T1_gauss * y.col(0);
                                for (size_t ii = 0; ii < nx; ii++)
                                    for (size_t p = 0; p < nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                            for (size_t q = 0; q < ny; q++)
                                            {
                                                double t2 = pow(-1, q);
                                                for (size_t j = 1; j < ny-1; j++)
                                                    for (size_t i = 1; i < nx-1; i++)
                                                    {
                                                        size_t k = i + j*nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        A(k, p+q*nx) -= gl_x[ii].weight*gl_y[jj].weight*dt1*t2*r(1)/pow(1 + x2_gl_w(jj), 2)/r3; 
                                                    }
                                            }
                                        }
                                    }
                            }
                            else if (c == 1) // East
                            {
                                arma::vec xw = T2_gauss * x.row(nx-1).t();
                                arma::vec yw =-T2_gauss * y.row(nx-1).t();
                                for (size_t jj = 0; jj < ny; jj++)
                                    for (size_t q = 0; q < ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double I = gl_x[ii].weight*gl_y[jj].weight*dt2*r(1)/pow(1 - x1_gl_w(ii), 2)/r3;
                                                    for (size_t p = 0; p < nx; p++)
                                                        A(k, p+q*nx) -= I;
                                                }
                                        }
                                    }
                            }
                            else if (c == 2) // North
                            {
                                arma::vec xw = T1_gauss * x.col(ny-1);
                                arma::vec yw =-T1_gauss * y.col(ny-1);
                                for (size_t ii = 0; ii < nx; ii++)
                                    for (size_t p = 0; p < nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double I =-gl_x[ii].weight*gl_y[jj].weight*dt1*r(1)/pow(1 - x2_gl_w(jj), 2)/r3;
                                                    for (size_t q = 0; q < ny; q++)
                                                        A(k, p+q*nx) -= I;
                                                }
                                        }
                                    }
                            }
                            else // West
                            {
                                arma::vec xw = T2_gauss * x.row(0).t();
                                arma::vec yw =-T2_gauss * y.row(0).t();
                                for (size_t jj = 0; jj < ny; jj++)
                                    for (size_t q = 0; q < ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                            for (size_t p = 0; p < nx; p++)
                                            {
                                                double t1 = pow(-1, p);
                                                for (size_t j = 1; j < ny-1; j++)
                                                    for (size_t i = 1; i < nx-1; i++)
                                                    {
                                                        size_t k = i + j*nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        A(k, p+q*nx) += gl_x[ii].weight*gl_y[jj].weight*t1*dt2*r(1)/pow(1 + x1_gl_w(ii), 2)/r3;
                                                    }
                                            }
                                        }
                                    }
                            }
                        }
            }
            break;
        }
        case Analysis::nonlinear:
        {
            arma::mat dzdxi_1 = D1*zC;
            arma::mat dzdxi_2 = zC*D2.t();
            arma::mat d2zdxi_12 = D1*dzdxi_1;
            arma::mat d2zdxi_1dxi_2 = D1*dzdxi_2;
            arma::mat d2zdxi_22 = dzdxi_2*D2.t();

            arma::mat dnxdxi_1 = D1*reshape(nC.col(0), nx, ny);
            arma::mat dnxdxi_2 = reshape(nC.col(0), nx, ny)*D2.t();
            arma::mat dnydxi_1 = D1*reshape(nC.col(1), nx, ny);
            arma::mat dnydxi_2 = reshape(nC.col(1), nx, ny)*D2.t();
            arma::mat dnzdxi_1 = D1*reshape(nC.col(2), nx, ny);
            arma::mat dnzdxi_2 = reshape(nC.col(2), nx, ny)*D2.t();
            #pragma omp parallel for reduction(+:A)
            for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
            {
                double y_lower = std::max(-1., xi_2(j)-delta/2);
                double y_upper = std::min(xi_2(j)+delta/2, 1.);
                for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                {
                    double x_left  = std::max(-1., xi_1(i)-delta/2);
                    double x_right = std::min(xi_1(i)+delta/2, 1.);

                    arma::vec::fixed<3> dXdxi_1       = {J11(i, j), J21(i, j), dzdxi_1(i, j)};
                    arma::vec::fixed<3> dXdxi_2       = {J12(i, j), J22(i, j), dzdxi_2(i, j)};
                    arma::vec::fixed<3> d2Xdxi_12     = {dJ11dxi_1(i, j), dJ21dxi_1(i, j), d2zdxi_12(i, j)};
                    arma::vec::fixed<3> d2Xdxi_1dxi_2 = {dJ11dxi_2(i, j), dJ22dxi_1(i, j), d2zdxi_1dxi_2(i, j)};
                    arma::vec::fixed<3> d2Xdxi_22     = {dJ12dxi_2(i, j), dJ22dxi_2(i, j), d2zdxi_22(i, j)};

                    arma::vec rho_tilde = externalContour(x_left, x_right, y_lower, y_upper, xi_1(i), xi_2(j), theta); 

                    size_t k = i+j*nx;
                    double d2T2   = 0;
                    double d2T2p1 = 0;
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double d2T1   = 0;
                        double d2T1p1 = 0;
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double m1  = dT1(i, p) *  T2(j, q);
                            double m2  =  T1(i, p) * dT2(j, q);
                            double m11 =      d2T1 *  T2(j, q);
                            double m12 = dT1(i, p) * dT2(j, q);
                            double m22 =  T1(i, p) * d2T2;
                            double dmudx = (J22(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_2(i, j))*m1
                                         - (J21(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_1(i, j))*m2;
                            double dmudy = (J11(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_1(i, j))*m2
                                         - (J12(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_2(i, j))*m1;
                            double dmudz = (J12(i, j)*nC(k, 1) - nC(k, 0)*J22(i, j))*m1
                                         - (J11(i, j)*nC(k, 1) - nC(k, 0)*J21(i, j))*m2;
                            arma::vec::fixed<3> gamma_0 = {dmudy*nC(k, 2) - dmudz*nC(k, 1),
                                                           dmudz*nC(k, 0) - dmudx*nC(k, 2),
                                                           dmudx*nC(k, 1) - dmudy*nC(k, 0)};
                            double dmudxdxi_1 = (J22(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_2(i, j))*m11 - (J21(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_1(i, j))*m12
                                              + (dJ22dxi_1(i, j)*nC(k, 2) + J22(i, j)*dnzdxi_1(i, j) - dnydxi_1(i, j)*dzdxi_2(i, j) - nC(k, 1)*d2zdxi_1dxi_2(i, j))*m1
                                              - (dJ21dxi_1(i, j)*nC(k, 2) + J21(i, j)*dnzdxi_1(i, j) - dnydxi_1(i, j)*dzdxi_1(i, j) - nC(k, 1)*d2zdxi_12(i, j))*m2;
                            double dmudxdxi_2 = (J22(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_2(i, j))*m12 - (J21(i, j)*nC(k, 2) - nC(k, 1)*dzdxi_1(i, j))*m22
                                              + (dJ22dxi_2(i, j)*nC(k, 2) + J22(i, j)*dnzdxi_2(i, j) - dnydxi_2(i, j)*dzdxi_2(i, j) - nC(k, 1)*d2zdxi_22(i, j))*m1
                                              - (dJ22dxi_1(i, j)*nC(k, 2) + J21(i, j)*dnzdxi_2(i, j) - dnydxi_2(i, j)*dzdxi_1(i, j) - nC(k, 1)*d2zdxi_1dxi_2(i, j))*m2;
                            double dmudydxi_1 = (J11(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_1(i, j))*m12 - (J12(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_2(i, j))*m11
                                              + (dJ11dxi_1(i, j)*nC(k, 2) + J11(i, j)*dnzdxi_1(i, j) - dnxdxi_1(i, j)*dzdxi_1(i, j) - nC(k, 0)*d2zdxi_12(i, j))*m2
                                              - (dJ11dxi_2(i, j)*nC(k, 2) + J12(i, j)*dnzdxi_1(i, j) - dnxdxi_1(i, j)*dzdxi_2(i, j) - nC(k, 0)*d2zdxi_1dxi_2(i, j))*m1;
                            double dmudydxi_2 = (J11(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_1(i, j))*m22 - (J12(i, j)*nC(k, 2) - nC(k, 0)*dzdxi_2(i, j))*m12
                                              + (dJ11dxi_2(i, j)*nC(k, 2) + J11(i, j)*dnzdxi_2(i, j) - dnxdxi_2(i, j)*dzdxi_1(i, j) - nC(k, 0)*d2zdxi_1dxi_2(i, j))*m2
                                              - (dJ12dxi_2(i, j)*nC(k, 2) + J12(i, j)*dnzdxi_2(i, j) - dnxdxi_2(i, j)*dzdxi_2(i, j) - nC(k, 0)*d2zdxi_22(i, j))*m1;
                            double dmudzdxi_1 = (J12(i, j)*nC(k, 1) - nC(k, 0)*J22(i, j))*m11 - (J11(i, j)*nC(k, 1) - nC(k, 0)*J21(i, j))*m12
                                              + (dJ11dxi_2(i, j)*nC(k, 1) + J12(i, j)*dnydxi_1(i, j) - dnxdxi_1(i, j)*J22(i, j) - nC(k, 0)*dJ22dxi_1(i, j))*m1
                                              - (dJ11dxi_1(i, j)*nC(k, 1) + J11(i, j)*dnydxi_1(i, j) - dnxdxi_1(i, j)*J21(i, j) - nC(k, 0)*dJ21dxi_1(i, j))*m2;
                            double dmudzdxi_2 = (J12(i, j)*nC(k, 1) - nC(k, 0)*J22(i, j))*m12 - (J11(i, j)*nC(k, 1) - nC(k, 0)*J21(i, j))*m22
                                              + (dJ12dxi_2(i, j)*nC(k, 1) + J12(i, j)*dnydxi_2(i, j) - dnxdxi_2(i, j)*J22(i, j) - nC(k, 0)*dJ22dxi_2(i, j))*m1
                                              - (dJ11dxi_2(i, j)*nC(k, 1) + J11(i, j)*dnydxi_2(i, j) - dnxdxi_2(i, j)*J21(i, j) - nC(k, 0)*dJ22dxi_1(i, j))*m2;
                            arma::vec::fixed<3> gamma_1_1 = {dmudydxi_1*nC(k, 2) + dmudy*dnzdxi_1(i, j) - dmudzdxi_1*nC(k, 1) - dmudz*dnydxi_1(i, j),
                                                             dmudzdxi_1*nC(k, 0) + dmudz*dnxdxi_1(i, j) - dmudxdxi_1*nC(k, 2) - dmudx*dnzdxi_1(i, j),
                                                             dmudxdxi_1*nC(k, 1) + dmudx*dnydxi_1(i, j) - dmudydxi_1*nC(k, 0) - dmudy*dnxdxi_1(i, j)};
                            arma::vec::fixed<3> gamma_1_2 = {dmudydxi_2*nC(k, 2) + dmudy*dnzdxi_2(i, j) - dmudzdxi_2*nC(k, 1) - dmudz*dnydxi_2(i, j),
                                                             dmudzdxi_2*nC(k, 0) + dmudz*dnxdxi_2(i, j) - dmudxdxi_2*nC(k, 2) - dmudx*dnzdxi_2(i, j),
                                                             dmudxdxi_2*nC(k, 1) + dmudx*dnydxi_2(i, j) - dmudydxi_2*nC(k, 0) - dmudy*dnxdxi_2(i, j)};
                            arma::mat gamma_1 = gamma_1_1*cT.t() + gamma_1_2*sT.t();
                            arma::vec::fixed<3> I0(arma::fill::zeros), I_1(arma::fill::zeros);
                            for (size_t n = 0; n < n_theta-1; n++)
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
                                I0  += F0*rho_tilde(n) + F1*pow(rho_tilde(n), 2)/2 + F2*pow(rho_tilde(n), 3)/3;
                                I_1 += F_1*log(fabs(rho_tilde(n)*A_abs));
                            }
                            arma::vec::fixed<3> q_mu = dtheta*(I0 + I_1);
                            A(k, p+q*nx) = dot(q_mu, nC.row(k));
                            std::swap(d2T1, d2T1p1);
                            if (p == 0)
                                d2T1p1 = 4;
                            else if (p < nx-1)
                                d2T1p1 = (p+2)*(2*dT1(i, p+1) + d2T1p1/p);
                        }
                        std::swap(d2T2, d2T2p1);
                        if (q == 0)
                            d2T2p1 = 4;
                        else if (q < ny-1)
                            d2T2p1 = (q+2)*(2*dT2(j, q+1) + d2T2p1/q);
                    }
                    if (y_lower > -1)
                        regularIntegralNonlinear(k, xC(i, j), yC(i, j), zC(i, j), std::ceil((x_right+1)/2*nx), std::ceil((y_lower+1)/2*ny), -1, x_right, -1, y_lower);
                    if (y_upper < 1)
                        regularIntegralNonlinear(k, xC(i, j), yC(i, j), zC(i, j), std::ceil((1-x_left )/2*nx), std::ceil((1-y_upper)/2*ny),  x_left,  1, y_upper,  1);
                    if (x_left > -1)
                        regularIntegralNonlinear(k, xC(i, j), yC(i, j), zC(i, j), std::ceil((x_left+1 )/2*nx), std::ceil((1-y_lower)/2*ny), -1,  x_left, y_lower,  1);
                    if (x_right < 1)
                        regularIntegralNonlinear(k, xC(i, j), yC(i, j), zC(i, j), std::ceil((1-x_right)/2*nx), std::ceil((y_upper+1)/2*ny), x_right,  1, -1, y_upper);
                }
            }
            for (const Wake* w:wakes)
                for (size_t c = 0; c < 4; c++)
                    if (chi[c] == w->chi)
                    {
                        std::vector<fastgl::QuadPair> gl_x_w(nx), gl_y_w(ny);
                        arma::vec x1_gl_w(nx, arma::fill::none), x2_gl_w(ny, arma::fill::none);
                        for (size_t ii = 0; ii < nx; ii++)
                        {
                            gl_x_w[ii] = fastgl::GLPair(nx, ii+1);
                            x1_gl_w(ii) =-gl_x_w[ii].x();
                        }
                        for (size_t jj = 0; jj < ny; jj++)
                        {
                            gl_y_w[jj] = fastgl::GLPair(ny, jj+1);
                            x2_gl_w(jj) =-gl_y_w[jj].x();
                        }
                        arma::vec x1_gauss = Chebyshev::gaussLobatto(nx);
                        arma::vec x2_gauss = Chebyshev::gaussLobatto(ny);
                        arma::mat D1_gauss = Chebyshev::derivativeMatrix(x1_gauss, Derivative::first);
                        arma::mat D2_gauss = Chebyshev::derivativeMatrix(x2_gauss, Derivative::first);
                        arma::mat dxdx1_gauss = D1_gauss * x;
                        arma::mat dydx1_gauss = D1_gauss * y;
                        arma::mat dzdx1_gauss = D1_gauss * z;
                        arma::mat dxdx2_gauss = x * D2_gauss.t();
                        arma::mat dydx2_gauss = y * D2_gauss.t();
                        arma::mat dzdx2_gauss = z * D2_gauss.t();
                        arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                        arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                        if (c == 0) // South
                        {
                            arma::vec xw    = T1_gauss *           x.col(0);
                            arma::vec yw    = T1_gauss *           y.col(0);
                            arma::vec zw    = T1_gauss *           z.col(0);
                            arma::vec dxdx1 = T1_gauss * dxdx1_gauss.col(0);
                            arma::vec dydx1 = T1_gauss * dydx1_gauss.col(0);
                            arma::vec dzdx1 = T1_gauss * dzdx1_gauss.col(0);
                            for (size_t ii = 0; ii < nx; ii++)
                            {
                                arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                for (size_t p = 0; p < nx; p++)
                                {
                                    double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                    for (size_t jj = 0; jj < ny; jj++)
                                    {
                                        double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                        double zW = zw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))*tan(alpha)/2;
                                        double dxdx2 =-pow(1 + x2_gl_w(jj),-2);
                                        // dydx2 = 0;
                                        double dzdx2 =-tan(alpha)/pow(1 + x2_gl_w(jj), 2);
                                        arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                        double e_11 = dot(e_1, e_1);
                                        double e_12 = dot(e_1, e_2);
                                        double e_22 = dot(e_2, e_2);
                                        double e = e_11*e_22 - pow(e_12, 2);
                                        double e11 = e_22/e;
                                        double e12 =-e_12/e;
                                        double e22 = e_11/e;
                                        for (size_t q = 0; q < ny; q++)
                                        {
                                            double t2 = pow(-1, q);
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = dt1 * t2;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                   dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                  -dydx1(ii)*dxdx2})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {-n_W(1)*dzdx2,
                                                                                  n_W(0)*dzdx2 - dxdx2*n_W(2),
                                                                                  dxdx2*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    A(k, p+q*nx) += dot(q_mu, nC.row(k)); 
                                                }
                                        }
                                    }
                                }
                            }
                        }
                        else if (c == 1) // East
                        {
                            arma::vec xw    = T2_gauss *           x.row(nx-1).t();
                            arma::vec yw    = T2_gauss *           y.row(nx-1).t();
                            arma::vec zw    = T2_gauss *           z.row(nx-1).t();
                            arma::vec dxdx2 = T2_gauss * dxdx2_gauss.row(nx-1).t();
                            arma::vec dydx2 = T2_gauss * dydx2_gauss.row(nx-1).t();
                            arma::vec dzdx2 = T2_gauss * dzdx2_gauss.row(nx-1).t();
                            for (size_t jj = 0; jj < ny; jj++)
                            {
                                arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                for (size_t q = 0; q < ny; q++)
                                {
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                    for (size_t ii = 0; ii < nx; ii++)
                                    {
                                        double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                        double zW = zw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))*tan(alpha)/2;
                                        double dxdx1 = 1/pow(1 - x1_gl_w(ii), 2);
                                        // dydx1 = 0
                                        double dzdx1 = tan(alpha)/pow(1 - x1_gl_w(ii), 2);
                                        arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                        double e_11 = dot(e_1, e_1);
                                        double e_12 = dot(e_1, e_2);
                                        double e_22 = dot(e_2, e_2);
                                        double e = e_11*e_22 - pow(e_12, 2);
                                        double e11 = e_22/e;
                                        double e12 =-e_12/e;
                                        double e22 = e_11/e;
                                        for (size_t j = 1; j < ny-1; j++)
                                            for (size_t i = 1; i < nx-1; i++)
                                            {
                                                size_t k = i + j*nx;
                                                arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                double r3 = pow(norm(r), 3);
                                                double dmudxi = dt2;
                                                double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                dxdx1*dydx2(jj)})/sqrt_A;
                                                arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                             dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                            -dxdx1*n_W(1)};
                                                arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                for (size_t p = 0; p < nx; p++)
                                                    A(k, p+q*nx) += dot(q_mu, nC.row(k));
                                            }
                                    }
                                }
                            }
                        }
                        else if (c == 2) // North
                        {
                            arma::vec xw    = T1_gauss *           x.col(ny-1);
                            arma::vec yw    = T1_gauss *           y.col(ny-1);
                            arma::vec zw    = T1_gauss *           z.col(ny-1);
                            arma::vec dxdx1 = T1_gauss * dxdx1_gauss.col(ny-1);
                            arma::vec dydx1 = T1_gauss * dydx1_gauss.col(ny-1);
                            arma::vec dzdx1 = T1_gauss * dzdx1_gauss.col(ny-1);
                            for (size_t ii = 0; ii < nx; ii++)
                            {
                                arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                for (size_t p = 0; p < nx; p++)
                                {
                                    double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                    for (size_t jj = 0; jj < ny; jj++)
                                    {
                                        double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                        double zW = zw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))*tan(alpha)/2;
                                        double dxdx2 = 1/pow(1 - x2_gl_w(jj), 2);
                                        // dydx2 = 0
                                        double dzdx2 = tan(alpha)/pow(1 - x2_gl_w(jj), 2);
                                        arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                        double e_11 = dot(e_1, e_1);
                                        double e_12 = dot(e_1, e_2);
                                        double e_22 = dot(e_2, e_2);
                                        double e = e_11*e_22 - pow(e_12, 2);
                                        double e11 = e_22/e;
                                        double e12 =-e_12/e;
                                        double e22 = e_11/e;
                                        for (size_t j = 1; j < ny-1; j++)
                                            for (size_t i = 1; i < nx-1; i++)
                                            {
                                                size_t k = i + j*nx;
                                                arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                double r3 = pow(norm(r), 3);
                                                double dmudxi = dt1;
                                                double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                               dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                              -dydx1(ii)*dxdx2})/sqrt_A;
                                                arma::vec::fixed<3> J_red = {-dzdx2*n_W(1),
                                                                              dzdx2*n_W(0) - dxdx2*n_W(2),
                                                                              dxdx2*n_W(1)};
                                                arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                for (size_t q = 0; q < ny; q++)
                                                    A(k, p+q*nx) += dot(q_mu, nC.row(k));
                                            }
                                    }
                                }
                            }
                        }
                        else // West
                        {
                            arma::vec xw    = T2_gauss *           x.row(0).t();
                            arma::vec yw    = T2_gauss *           y.row(0).t();
                            arma::vec zw    = T2_gauss *           z.row(0).t();
                            arma::vec dxdx2 = T2_gauss * dxdx2_gauss.row(0).t();
                            arma::vec dydx2 = T2_gauss * dydx2_gauss.row(0).t();
                            arma::vec dzdx2 = T2_gauss * dzdx2_gauss.row(0).t();
                            for (size_t jj = 0; jj < ny; jj++)
                            {
                                arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                for (size_t q = 0; q < ny; q++)
                                {
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                    for (size_t ii = 0; ii < nx; ii++)
                                    {
                                        double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                        double zW = zw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))*tan(alpha)/2;
                                        double dxdx1 =-1/pow(1 + x1_gl_w(ii), 2);
                                        // dydx1 = 0
                                        double dzdx1 =-tan(alpha)/pow(1 + x1_gl_w(ii), 2);
                                        arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                        double e_11 = dot(e_1, e_1);
                                        double e_12 = dot(e_1, e_2);
                                        double e_22 = dot(e_2, e_2);
                                        double e = e_11*e_22 - pow(e_12, 2);
                                        double e11 = e_22/e;
                                        double e12 =-e_12/e;
                                        double e22 = e_11/e;
                                        for (size_t p = 0; p < nx; p++)
                                        {
                                            double t1 = pow(-1, p);
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = t1*dt2;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                    dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                    dxdx1*dydx2(jj)})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                 dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                -dxdx1*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    A(k, p+q*nx) += dot(q_mu, nC.row(k)); 
                                                }
                                        }
                                    }
                                }
                            }
                        }
                    }
            if (sym == Symmetry::y)
            {
                std::vector<fastgl::QuadPair> gl_x(nx), gl_y(ny);
                arma::vec x1_gl(nx, arma::fill::none), x2_gl(ny, arma::fill::none);
                for (size_t ii = 0; ii < nx; ii++)
                {
                    gl_x[ii] = fastgl::GLPair(nx, ii+1);
                    x1_gl(ii) =-gl_x[ii].x();
                }
                for (size_t jj = 0; jj < ny; jj++)
                {
                    gl_y[jj] = fastgl::GLPair(ny, jj+1);
                    x2_gl(jj) =-gl_y[jj].x();
                }
                arma::mat Tx_gl = Lagrange::interpolationMatrix(x1, x1_gl);
                arma::mat Ty_gl = Lagrange::interpolationMatrix(x2, x2_gl);
                arma::mat z_gl = Lagrange::interpolation2D(Tx_gl, Ty_gl, z, x1_gl, x2_gl);
                arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat dz_gldx1 = D1_gl*z_gl;
                arma::mat dz_gldx2 = z_gl*D2_gl.t();
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
                y_gl     =-y_gl;
                dy_gldx1 =-dy_gldx1;
                dy_gldx2 =-dy_gldx2;
                arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
                arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
                #pragma omp parallel for reduction(+:A)
                for (size_t j = 1; j < ny-1; j++) // Loop over Collocation Points in 2-direction
                    for (size_t i = 1; i < nx-1; i++) // Loop over Collocation Points in 1-direction
                        for (size_t ii = 0; ii < nx; ii++)
                            for (size_t jj = 0; jj < ny; jj++)
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
                                double r3 = pow(norm(r), 3);
                                arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                                                dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                                                dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
                                arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                                                {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                                                {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
                                double t2   = 1; // First Chebyshev Polynomial (j)
                                double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

                                double dt2   = 0;
                                double dt2p1 = 1;
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                {
                                    double t1   = 1; // First Chebyshev Polynomial (i)
                                    double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                                    double dt1   = 0;
                                    double dt1p1 = 1;
                                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                    {
                                        double T_1 = dt1 * t2;
                                        double T_2 =  t1 * dt2;
                                        arma::vec::fixed<2> dmudxi = {T_1, T_2};
                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                        arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3;
                                        A(i+j*nx, p+q*nx) -= dot(q_mu, nC.row(i+j*nx));
                                        std::swap(t1, t1p1);
                                        t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1);

                                        std::swap(dt1, dt1p1);
                                        dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                                    }
                                    std::swap(t2, t2p1);
                                    t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1);

                                    std::swap(dt2, dt2p1);
                                    dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
                                }
                            }
                for (const Wake* w:wakes)
                    for (size_t c = 0; c < 4; c++)
                        if (chi[c] == w->chi)
                        {
                            std::vector<fastgl::QuadPair> gl_x_w(nx), gl_y_w(ny);
                            arma::vec x1_gl_w(nx, arma::fill::none), x2_gl_w(ny, arma::fill::none);
                            for (size_t ii = 0; ii < nx; ii++)
                            {
                                gl_x_w[ii] = fastgl::GLPair(nx, ii+1);
                                x1_gl_w(ii) =-gl_x_w[ii].x();
                            }
                            for (size_t jj = 0; jj < ny; jj++)
                            {
                                gl_y_w[jj] = fastgl::GLPair(ny, jj+1);
                                x2_gl_w(jj) =-gl_y_w[jj].x();
                            }
                            arma::vec x1_gauss = Chebyshev::gauss(nx);
                            arma::vec x2_gauss = Chebyshev::gauss(ny);
                            arma::mat D1_gauss = Lagrange::derivativeMatrix(x1_gauss);
                            arma::mat D2_gauss = Lagrange::derivativeMatrix(x2_gauss);
                            arma::mat dxdx1_gauss = D1_gauss * x;
                            arma::mat dydx1_gauss = D1_gauss * y;
                            arma::mat dzdx1_gauss = D1_gauss * z;
                            arma::mat dxdx2_gauss = x * D2_gauss.t();
                            arma::mat dydx2_gauss = y * D2_gauss.t();
                            arma::mat dzdx2_gauss = z * D2_gauss.t();
                            arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                            arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                            if (c == 0) // South
                            {
                                arma::vec xw    = T1_gauss *           x.col(0);
                                arma::vec yw    =-T1_gauss *           y.col(0);
                                arma::vec zw    = T1_gauss *           z.col(0);
                                arma::vec dxdx1 = T1_gauss * dxdx1_gauss.col(0);
                                arma::vec dydx1 =-T1_gauss * dydx1_gauss.col(0);
                                arma::vec dzdx1 = T1_gauss * dzdx1_gauss.col(0);
                                for (size_t ii = 0; ii < nx; ii++)
                                {
                                    arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                    for (size_t p = 0; p < nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                            double zW = zw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))*tan(alpha)/2;
                                            double dxdx2 =-pow(1 + x2_gl_w(jj),-2);
                                            // dydx2 = 0;
                                            double dzdx2 =-tan(alpha)/pow(1 + x2_gl_w(jj), 2);
                                            arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t q = 0; q < ny; q++)
                                            {
                                                double t2 = pow(-1, q);
                                                for (size_t j = 1; j < ny-1; j++)
                                                    for (size_t i = 1; i < nx-1; i++)
                                                    {
                                                        size_t k = i + j*nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = dt1 * t2;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                       dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                      -dydx1(ii)*dxdx2})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {-n_W(1)*dzdx2,
                                                                                        n_W(0)*dzdx2 - dxdx2*n_W(2),
                                                                                        dxdx2*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        A(k, p+q*nx) -= dot(q_mu, nC.row(k)); 
                                                    }
                                            }
                                        }
                                    }
                                }
                            }
                            else if (c == 1) // East
                            {
                                arma::vec xw    = T2_gauss *           x.row(nx-1).t();
                                arma::vec yw    =-T2_gauss *           y.row(nx-1).t();
                                arma::vec zw    = T2_gauss *           z.row(nx-1).t();
                                arma::vec dxdx2 = T2_gauss * dxdx2_gauss.row(nx-1).t();
                                arma::vec dydx2 =-T2_gauss * dydx2_gauss.row(nx-1).t();
                                arma::vec dzdx2 = T2_gauss * dzdx2_gauss.row(nx-1).t();
                                for (size_t jj = 0; jj < ny; jj++)
                                {
                                    arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                    for (size_t q = 0; q < ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                            double zW = zw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))*tan(alpha)/2;
                                            double dxdx1 = 1/pow(1 - x1_gl_w(ii), 2);
                                            // dydx1 = 0
                                            double dzdx1 = tan(alpha)/pow(1 - x1_gl_w(ii), 2);
                                            arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = dt2;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                    dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                    dxdx1*dydx2(jj)})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                 dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                -dxdx1*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    for (size_t p = 0; p < nx; p++)
                                                        A(k, p+q*nx) -= dot(q_mu, nC.row(k));
                                                }
                                        }
                                    }
                                }
                            }
                            else if (c == 2) // North
                            {
                                arma::vec xw    = T1_gauss *           x.col(ny-1);
                                arma::vec yw    =-T1_gauss *           y.col(ny-1);
                                arma::vec zw    = T1_gauss *           z.col(ny-1);
                                arma::vec dxdx1 = T1_gauss * dxdx1_gauss.col(ny-1);
                                arma::vec dydx1 =-T1_gauss * dydx1_gauss.col(ny-1);
                                arma::vec dzdx1 = T1_gauss * dzdx1_gauss.col(ny-1);
                                for (size_t ii = 0; ii < nx; ii++)
                                {
                                    arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                    for (size_t p = 0; p < nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                            double zW = zw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))*tan(alpha)/2;
                                            double dxdx2 = 1/pow(1 - x2_gl_w(jj), 2);
                                            // dydx2 = 0
                                            double dzdx2 = tan(alpha)/pow(1 - x2_gl_w(jj), 2);
                                            arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t j = 1; j < ny-1; j++)
                                                for (size_t i = 1; i < nx-1; i++)
                                                {
                                                    size_t k = i + j*nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = dt1;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                   dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                  -dydx1(ii)*dxdx2})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {-dzdx2*n_W(1),
                                                                                  dzdx2*n_W(0) - dxdx2*n_W(2),
                                                                                  dxdx2*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    for (size_t q = 0; q < ny; q++)
                                                        A(k, p+q*nx) -= dot(q_mu, nC.row(k));
                                                }
                                        }
                                    }
                                }
                            }
                            else // West
                            {
                                arma::vec xw    = T2_gauss *           x.row(0).t();
                                arma::vec yw    =-T2_gauss *           y.row(0).t();
                                arma::vec zw    = T2_gauss *           z.row(0).t();
                                arma::vec dxdx2 = T2_gauss * dxdx2_gauss.row(0).t();
                                arma::vec dydx2 =-T2_gauss * dydx2_gauss.row(0).t();
                                arma::vec dzdx2 = T2_gauss * dzdx2_gauss.row(0).t();
                                for (size_t jj = 0; jj < ny; jj++)
                                {
                                    arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                    for (size_t q = 0; q < ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                            double zW = zw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))*tan(alpha)/2;
                                            double dxdx1 =-1/pow(1 + x1_gl_w(ii), 2);
                                            // dydx1 = 0
                                            double dzdx1 =-tan(alpha)/pow(1 + x1_gl_w(ii), 2);
                                            arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t p = 0; p < nx; p++)
                                            {
                                                double t1 = pow(-1, p);
                                                for (size_t j = 1; j < ny-1; j++)
                                                    for (size_t i = 1; i < nx-1; i++)
                                                    {
                                                        size_t k = i + j*nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = t1*dt2;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                        dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                        dxdx1*dydx2(jj)})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                        dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                        -dxdx1*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        A(k, p+q*nx) -= dot(q_mu, nC.row(k)); 
                                                    }
                                            }
                                        }
                                    }
                                }
                            }
                        }
            }
            break;
        }
        default:
            std::println("Only linear and nonlinear analysis are implemented for Wing!");
            exit(EXIT_FAILURE);
    }
}