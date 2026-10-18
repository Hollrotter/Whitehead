#include "Wing.hpp"

void Wing::integralRightLinear(size_t i, size_t j, size_t n1, size_t n2, double x_right)
{

}

void Wing::integralLeftLinear(size_t i, size_t j, size_t n1, size_t n2, double x_left)
{

}

void Wing::integralUpperLinear(size_t i, size_t j, size_t n1, size_t n2, double y_upper)
{

}

void Wing::integralLowerLinear(size_t i, size_t j, size_t n1, size_t n2, double y_lower)
{

}

void Wing::integralLowerRightLinear(size_t i, size_t j, size_t n1, size_t n2, double x_right, double y_upper)
{

}

void Wing::integralLowerLeftLinear(size_t i, size_t j, size_t n1, size_t n2, double x_right, double y_lower)
{

}

void Wing::integralUpperRightLinear(size_t i, size_t j, size_t n1, size_t n2, double x_left, double y_upper)
{

}

void Wing::integralUpperLeftLinear(size_t i, size_t j, size_t n1, size_t n2, double x_left, double y_lower)
{

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
            for (size_t ii = 0; ii < nx_ii; ii++)
            {
                arma::vec::fixed<2> r = {x_gl(ii, jj) - xC, y_gl(ii, jj) - yC};
                double r3 = pow(norm(r), 3);

                double  w1 = phi1->weightFunction(x1_gl(ii));
                double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));

                double t1   = phi1->constant(); // First Chebyshev Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Chebyshev Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    double dmudx1 = dpsi1 *  psi2;
                    double dmudx2 =  psi1 * dpsi2;
                    A(k, p+q*nx) += gl_x[ii].weight*gl_y[jj].weight/r3*(xp - xm)*(yp - ym)/4
                                    *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                    - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                    std::swap(t1, t1p1); // Swap order so the oldes Chebyshev Polynomial will be overwritten
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
            }
            std::swap(t2, t2p1); // Swap order so the oldest Chebyshev Polynomial will be overwritten
            phi2->next(q, x2_gl(jj), t2, t2p1);

            std::swap(dt2, dt2p1);
            phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
        }
    }
}