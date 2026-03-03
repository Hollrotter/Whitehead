#include "Splinefit.hpp"

/**
 * @brief 
 * 
 * @param x 
 */
void Splinefit::setQR(const arma::vec x)
{
    arma::mat S(x.size(), M-1+order);
    B = B_Spline(arma::linspace(x.front(), x.back(), M), order);
    #pragma omp parallel for
    for (size_t m = 0; m < M-1+order; m++)
        S.col(m) = B(x, m);
    arma::qr_econ(Q, R, S);
    Q = Q.t();
}