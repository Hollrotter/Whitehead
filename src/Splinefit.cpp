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
 * @brief Construct a new Splinefit:: Splinefit object
 * 
 * @param x 
 * @param y 
 * @param _M 
 * @param _N 
 */
Splinefit::Splinefit(arma::vec x, arma::vec y, size_t _M, size_t _N) : M(_M), N(_N)
{
    setQR(x, y);
}

/**
 * @brief Construct a new Splinefit:: Splinefit object
 * 
 * @param x 
 * @param y 
 * @param z 
 * @param _M 
 * @param _N 
 */
Splinefit::Splinefit(arma::vec x, arma::vec y, arma::mat z, size_t _M, size_t _N) : M(_M), N(_N)
{
    setQR(x, y);
    fit(z);
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
        y_dach += B1(x, j) * h(j);
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

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return double 
 */
double Splinefit::operator()(const double x, const double y)
{
    double z_dach = 0;
    #pragma omp parallel for
    for (size_t m = 0; m < M-1+order; m++)
        for (size_t n = 0; n < N-1+order; n++)
            z_dach += B1(x, m) * B2(y, n) * h(n+m*(N-1+order));
    return z_dach;
}

/**
 * @brief 
 * 
 * @param x 
 * @param y 
 * @return arma::mat 
 */
arma::mat Splinefit::operator()(const arma::vec x, const arma::vec y)
{
    arma::mat z_dach(x.size(), y.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = 0; j < y.size(); j++)
            z_dach(i, j) = operator()(x(i), y(j));
    return z_dach;
}