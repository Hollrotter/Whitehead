#include "Camber.hpp"

/**
 * @brief 
 * 
 * @param x Normalized x-coordinate where the camber height will be evaluated.
 * @return double 
 */
double Camber::operator() (double x)
{
    switch(camberType)
    {
        case CamberType::none:
            return 0;
        case CamberType::function:
            return f(x);
        case CamberType::b_spline:
            return s(x);
    }
    std::unreachable();
}

/**
 * @brief 
 * 
 * @param x Vector of normalized x-coordinates where the camber height will be evaluated.
 * @return arma::vec 
 */
arma::vec Camber::operator() (arma::vec x)
{
    arma::vec z(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        z(i) = operator()(x(i));
    return z;
}

/**
 * @brief 
 * 
 * @param x Normalized x-coordinate where the camber slope will be evaluated.
 * @return double 
 */
double Camber::diff(double x)
{
    switch(camberType)
    {
        case CamberType::none:
            return 0;
        case CamberType::function:
            return df(x);
        case CamberType::b_spline:
            return s.diff(x);
    }
    std::unreachable();
}

/**
 * @brief 
 * 
 * @param x Vector of normalized x-coordinates where the camber slope will be evaluated.
 * @return arma::vec 
 */
arma::vec Camber::diff(arma::vec x)
{
    arma::vec dz(x.size());
    #pragma omp parallel for
    for (size_t i = 0; i < x.size(); i++)
        dz(i) = diff(x(i));
    return dz;
}