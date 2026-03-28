#include "Airfoil.hpp"

/**
 * @brief 
 * 
 */
void Airfoil::aerodynamicMatrix()
{
        A.row(0).ones();
        for (size_t i = 1; i < nx; i++)
        {
            A(i, 0) = log((1 - xi(i))/(1 + xi(i)));
            A(i, 1) = 2 + xi(i)*A(i, 0);
            for (size_t j = 2; j < nx; j++)
                A(i, j) = 2*(Chebyshev::integral(j-1) + xi(i)*A(i, j-1)) - A(i, j-2);
        }
        if (analysis == Analysis::nonlinear)
        {
            nC.zeros(2, nx);
            auto [dxdxi, dzdxi] = chi->derivative(xi);
            arma::mat D = Lagrange::derivativeMatrix(xi);
            arma::vec d2xdxi2 = D*dxdxi;
            arma::vec d2zdxi2 = D*dzdxi;
            for (size_t i = 1; i < nx; i++)
            {
                double xi_min = i == 0    ? -1 : xi(i-1);
                double xi_max = i == nx-1 ?  1 : xi(i+1);
                double dr0 = sqrt(pow(dxdxi(i), 2) + pow(dzdxi(i), 2));
                double d2r0 = (dxdxi(i)*d2xdxi2(i) + dzdxi(i)*d2zdxi2(i))/dr0;
                nC.col(i) = arma::vec{-dzdxi(i), dxdxi(i)}/dr0;
                for (size_t j = 0; j < nx; j++)
                {
                    double F0 = boost::math::chebyshev_t(j, xi(i));
                    double F1 = boost::math::chebyshev_t_prime(j, xi(i));
                    double tau = (xi_max + xi_min)/2;
                    A(i, j) += d2r0/dr0*(F0 + F1*tau) * (xi_max - xi_min);
                    if (i > 0)
                    {
                        double I = 0;
                        for (size_t n = 0; n < nx-i; n++)
                        {
                            fastgl::QuadPair gl = fastgl::GLPair(nx-i, n+1);
                            double xi_gl = (1 - gl.x())/2*(xi_min + 1) - 1;
                            auto [ x_gl,  z_gl] = chi->evaluate(xi_gl);
                            auto [dx_gl, dz_gl] = chi->derivative(xi_gl);
                            double dr = sqrt(pow(dx_gl, 2) + pow(dz_gl, 2));
                            double T = boost::math::chebyshev_t(j, xi_gl);
                            I += gl.weight * k1(dr0, dr, dxdxi(i), dzdxi(i), x(i), x_gl, z(i), z_gl, xi(i), xi_gl) * T;
                        }
                        A(i, j) += I * (xi_min+1)/2;
                    }
                    if (i < nx-1)
                    {
                        double I = 0;
                        for (size_t n = 0; n < i+2; n++)
                        {
                            fastgl::QuadPair gl = fastgl::GLPair(i+2, n+1);
                            double xi_gl = (1 - gl.x())/2*(1 - xi_max) + xi_max;
                            auto [ x_gl,  z_gl] = chi->evaluate(xi_gl);
                            auto [dx_gl, dz_gl] = chi->derivative(xi_gl);
                            double dr = sqrt(pow(dx_gl, 2) + pow(dz_gl, 2));
                            double T = boost::math::chebyshev_t(j, xi_gl);
                            I += gl.weight * k1(dr0, dr, dxdxi(i), dzdxi(i), x(i), x_gl, z(i), z_gl, xi(i), xi_gl) * T;
                        }
                        A(i, j) += I * (1-xi_max)/2;
                    }
                }
            }
        }
}