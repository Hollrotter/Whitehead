#include "Wing.hpp"

void Wing::integralRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_right)
{

}

void Wing::integralLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_left)
{

}

void Wing::integralUpperNonlinear(size_t i, size_t j, size_t n1, size_t n2, double y_upper)
{

}

void Wing::integralLowerNonlinear(size_t i, size_t j, size_t n1, size_t n2, double y_lower)
{

}

void Wing::integralLowerRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_right, double y_upper)
{

}

void Wing::integralLowerLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_right, double y_lower)
{

}

void Wing::integralUpperRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_left, double y_upper)
{

}

void Wing::integralUpperLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double x_left, double y_lower)
{

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
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < ny_jj; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC, z_gl(ii, jj) - zC};
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
                    A(k, p+q*nx) += dot(q_mu, nC.row(k));
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