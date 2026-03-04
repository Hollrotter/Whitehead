#include "DVM.hpp"

/**
 * @brief 
 * 
 */
void DVM::aerodynamicMatrix()
{
    switch (analysis)
    {
        case Analysis::linear:
        {
            A = 1/arma::datum::tau/(repelem(xg.t(), nx, 1) - repelem(xC, 1, nx));
            break;
        }
        case Analysis::nonlinear:
        {
            #pragma omp parallel for
            for (size_t i = 0; i < nx; i++)
                for (size_t j = 0; j < nx; j++)
                    A(i, j) = dot(arma::vec::fixed<2>{zC(i)-zg(j), xg(j)-xC(i)}/arma::datum::tau/(pow(xg(j)-xC(i), 2) + pow(zg(j)-zC(i), 2)), nC.col(i));
            break;
        }
    }
}