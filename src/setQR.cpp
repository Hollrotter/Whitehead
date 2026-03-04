#include "Splinefit.hpp"

/**
 * @brief 
 * 
 * @param x 
 */
void Splinefit::setQR(const arma::vec x)
{
    arma::mat S(x.size(), M-1+order);
    B1 = B_Spline(arma::linspace(x.front(), x.back(), M), order);
    #pragma omp parallel for
    for (size_t m = 0; m < M-1+order; m++)
        S.col(m) = B1(x, m);
    arma::qr_econ(Q, R, S);
    Q = Q.t();
}

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 */
void Splinefit::setQR(const arma::vec x, const arma::vec y)
{
    arma::mat S(x.size()*y.size(), (M-1+order)*(N-1+order));
    B1 = B_Spline(arma::linspace(x.front(), x.back(), M), order);
    B2 = B_Spline(arma::linspace(y.front(), y.back(), N), order);
    #pragma omp parallel for
    for (size_t m = 0; m < M-1+order; m++)
        for (size_t n = 0; n < N-1+order; n++)
            S.col(n+m*(N-1+order)) = vectorise(B1(x, m)*B2(y, n).t());
    arma::qr_econ(Q, R, S);
    Q = Q.t();
}