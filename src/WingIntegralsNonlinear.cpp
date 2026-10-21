#include "Wing.hpp"

void Wing::integralRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp)
{
    arma::vec x1_gl(n1, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.5, 0.0);
        x1_gl = ((1-xp)*X1 + xp+1)/2;
        w1_gl = pow((1 - xp)/2, 1.5)*W1%sqrt(1+x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (1 - xp)/2 + xp;
            w1_gl(ii) = gl_x.weight*(1 - xp)/2;
        }
    arma::vec x2_gl = phi2->xg;
    arma::vec w2_gl = phi2->wg;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm)
{
    arma::vec x1_gl(n1, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.0, 0.5);
        x1_gl = ((xm+1)*X1 + xm-1)/2;
        w1_gl = pow((xm + 1)/2, 1.5)*W1%sqrt(1-x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (xm + 1)/2 - 1;
            w1_gl(ii) = gl_x.weight*(xm + 1)/2;
        }
    arma::vec x2_gl = phi2->xg;
    arma::vec w2_gl = phi2->wg;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralUpperNonlinear(size_t i, size_t j, size_t n1, size_t n2, double yp)
{
    arma::vec x1_gl = phi1->xg;
    arma::vec w1_gl = phi1->wg;
    arma::vec x2_gl(n2, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none);
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.5, 0.0);
        x2_gl = ((1-yp)*X2 + yp+1)/2;
        w2_gl = pow((1 - yp)/2, 1.5)*W2%sqrt(1+x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (1 - yp)/2 + yp;
            w2_gl(jj) = gl_y.weight*(1 - yp)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralLowerNonlinear(size_t i, size_t j, size_t n1, size_t n2, double ym)
{
    arma::vec x1_gl = phi1->xg;
    arma::vec w1_gl = phi1->wg;
    arma::vec x2_gl(n2, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none);
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.0, 0.5);
        x2_gl = ((ym+1)*X2 + ym-1)/2;
        w2_gl = pow((ym + 1)/2, 1.5)*W2%sqrt(1-x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (ym + 1)/2 - 1;
            w2_gl(jj) = gl_y.weight*(ym + 1)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralLowerRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp, double yp)
{
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), w2_gl(n2, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.5, 0.0);
        x1_gl = ((1-xp)*X1 + xp+1)/2;
        w1_gl = pow((1 - xp)/2, 1.5)*W1%sqrt(1+x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (1 - xp)/2 + xp;
            w1_gl(ii) = gl_x.weight*(1 - xp)/2;
        }
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.0, 0.5);
        x2_gl = ((yp+1)*X2 + yp-1)/2;
        w2_gl = pow((yp + 1)/2, 1.5)*W2%sqrt(1-x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (yp + 1)/2 - 1;
            w2_gl(jj) = gl_y.weight*(yp + 1)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralLowerLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp, double ym)
{
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), w2_gl(n2, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.0, 0.5);
        x1_gl = ((xp+1)*X1 + xp-1)/2;
        w1_gl = pow((xp + 1)/2, 1.5)*W1%sqrt(1-x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (xp + 1)/2 - 1;
            w1_gl(ii) = gl_x.weight*(xp + 1)/2;
        }
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.0, 0.5);
        x2_gl = ((ym+1)*X2 + ym-1)/2;
        w2_gl = pow((ym + 1)/2, 1.5)*W2%sqrt(1-x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (ym + 1)/2 - 1;
            w2_gl(jj) = gl_y.weight*(ym + 1)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralUpperRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double yp)
{
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), w2_gl(n2, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.5, 0.0);
        x1_gl = ((1-xm)*X1 + xm+1)/2;
        w1_gl = pow((1 - xm)/2, 1.5)*W1%sqrt(1+x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (1 - xm)/2 + xm;
            w1_gl(ii) = gl_x.weight*(1 - xm)/2;
        }
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.5, 0.0);
        x2_gl = ((1-yp)*X2 + yp+1)/2;
        w2_gl = pow((1 - yp)/2, 1.5)*W2%sqrt(1+x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (1 - yp)/2 + yp;
            w2_gl(jj) = gl_y.weight*(1 - yp)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::integralUpperLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double ym)
{
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), w2_gl(n2, arma::fill::none);
    if (phi1->basis == Basis::U)
    {
        auto [X1, W1] = gaujac(n1, 0.0, 0.5);
        x1_gl = ((xm+1)*X1 + xm-1)/2;
        w1_gl = pow((xm + 1)/2, 1.5)*W1%sqrt(1-x1_gl);
    }
    else
        for (size_t ii = 0; ii < n1; ii++)
        {
            fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
            x1_gl(ii) = (1 - gl_x.x()) * (xm + 1)/2 - 1;
            w1_gl(ii) = gl_x.weight*(xm + 1)/2;
        }
    if (phi2->basis == Basis::U)
    {
        auto [X2, W2] = gaujac(n2, 0.5, 0.0);
        x2_gl = ((1-ym)*X2 + ym+1)/2;
        w2_gl = pow((1 - ym)/2, 1.5)*W2%sqrt(1+x2_gl);
    }
    else
        for (size_t jj = 0; jj < n2; jj++)
        {
            fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
            x2_gl(jj) = (1 - gl_y.x()) * (1 - ym)/2 + ym;
            w2_gl(jj) = gl_y.weight*(1 - ym)/2;
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = w1_gl(ii) * w2_gl(jj) * cross(gamma_gl, r)/r3;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}

void Wing::regularIntegralNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double xp, double ym, double yp)
{
    std::vector<fastgl::QuadPair> gl_x(n1), gl_y(n2);
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    for (size_t ii = 0; ii < n1; ii++)
    {
        gl_x[ii] = fastgl::GLPair(n1, ii+1);
        x1_gl(ii) = (1 - gl_x[ii].x()) * (xp - xm)/2 + xm;
    }
    for (size_t jj = 0; jj < n2; jj++)
    {
        gl_y[jj] = fastgl::GLPair(n2, jj+1);
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
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Chebyshev Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Chebyshev Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3 * (xp - xm) * (yp - ym)/4;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1);
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1);
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}