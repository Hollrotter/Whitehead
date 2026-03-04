#include "Splinefit.hpp"

/**
 * @brief 
 * 
 * @param x 
 * @return double 
 */
double Splinefit::diff(const double x)
{
    double dydx = 0;
    #pragma omp parallel for
    for (size_t j = 0; j < h.size(); j++)
        dydx += B1.diff(x, j) * h(j);
    return dydx;
}

/**
 * @brief 
 * 
 * @param x 
 * @return arma::vec 
 */
arma::vec Splinefit::diff(const arma::vec x)
{
    arma::vec dydx(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        dydx(i) = diff(x(i));
    return dydx;
}

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return arma::vec 
 */
arma::vec Splinefit::diff(const double x, const double y)
{
    arma::vec gradZ(2);
    #pragma omp parallel for
    for (size_t m = 0; m < M-1+order; m++)
        for (size_t n = 0; n < N-1+order; n++)
            gradZ += arma::vec{B1.diff(x, m)*B2(y, n)*h(n+m*(N-1+order)), B1(x, m)*B2.diff(y, n)*h(n+m*(N-1+order))};
    return gradZ;
}

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return arma::cube 
 */
arma::cube Splinefit::diff(const arma::vec x, const arma::vec y)
{
    arma::cube gradZ(x.size(), y.size(), 2);
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = 0; j < y.size(); j++)
            gradZ.tube(i, j) = diff(x(i), y(j));
    return gradZ;
}