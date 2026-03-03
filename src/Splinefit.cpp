#include "Splinefit.hpp"

/**
 * @brief Construct a new Splinefit:: Splinefit object
 * 
 * @param x 
 * @param _M 
 */
Splinefit::Splinefit(arma::vec x, size_t _M) : M(_M)
{
    setQR(x);
}

/**
 * @brief Construct a new Splinefit:: Splinefit object
 * 
 * @param x 
 * @param y 
 * @param _M 
 */
Splinefit::Splinefit(arma::vec x, arma::vec y, size_t _M) : M(_M)
{
    setQR(x);
    fit(y);
}

/**
 * @brief 
 * 
 * @param x 
 * @return double 
 */
double Splinefit::operator()(const double x)
{
    double y_dach = 0;
    #pragma omp parallel for
    for (size_t j = 0; j < h.size(); j++)
        y_dach += B(x, j) * h(j);
    return y_dach;
}

/**
 * @brief 
 * 
 * @param x 
 * @return arma::vec 
 */
arma::vec Splinefit::operator()(const arma::vec x)
{
    arma::vec y_dach(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        y_dach(i) = operator()(x(i));
    return y_dach;
}