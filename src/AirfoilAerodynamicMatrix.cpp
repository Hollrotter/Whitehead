#include "Airfoil.hpp"

/**
 * @brief 
 * 
 */
void Airfoil::aerodynamicMatrix()
{
        for (size_t i = 0; i < nx; i++)
        {
            double v = 1; // V_j(xi): Chebyshev Polynomial of the third kind
            double vp1 = 2*xi(i) - 1; // V_j+1(xi): Next Chebyshev Polynomial
            for (size_t j = 0; j < nx; j++)
            {
                A(i, j) = v;
                std::swap(v, vp1);
                vp1 = boost::math::chebyshev_next(xi(i), v, vp1);
            }
        }
        if (analysis == Analysis::nonlinear)
        {
            nC.zeros(2, nx);
            auto [dxdxi, dzdxi] = chi->derivative(xi);
            arma::mat D = Lagrange::derivativeMatrix(xi);
            arma::vec d2xdxi2 = D*dxdxi;
            arma::vec d2zdxi2 = D*dzdxi;
            nC.col(0) = arma::vec::fixed<2>{-dzdxi(0), dxdxi(0)}/sqrt(pow(dxdxi(0), 2) + pow(dzdxi(0), 2));
            for (size_t i = 0; i < nx; i++)
            {
                double dr0 = sqrt(pow(dxdxi(i), 2) + pow(dzdxi(i), 2));
                nC.col(i) = arma::vec::fixed<2>{-dzdxi(i), dxdxi(i)}/dr0;
                for (size_t n = 0; n < nx; n++)
                {
                    double xi_w = cos(arma::datum::pi*(nx - n)/(nx+0.5));
                    auto [ x_gl,  z_gl] = chi->evaluate(xi_w);
                    auto [dx_gl, dz_gl] = chi->derivative(xi_w);
                    double dr = sqrt(pow(dx_gl, 2) + pow(dz_gl, 2));
                    double w = 1;
                    double wp1 = 2*xi_w + 1;
                    for (size_t j = 0; j < nx; j++)
                    {
                        A(i, j) += (1-xi_w)/(nx+0.5) * k1(dr0, dr, dxdxi(i), dzdxi(i), x(i), x_gl, z(i), z_gl, xi(i), xi_w) * w;
                        std::swap(w, wp1);
                        wp1 = boost::math::chebyshev_next(xi_w, w, wp1);
                    }
                }
            }
        }
}