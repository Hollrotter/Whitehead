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
        dydx += B.diff(x, j) * h(j);
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