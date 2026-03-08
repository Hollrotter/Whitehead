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
            break;
        }
    }
}