#include "misc.hpp"

/**
 * @brief 
 * 
 * @param a First value for comparison.
 * @param b Second value for comparison.
 * @return true 
 * @return false 
 */
bool almostEqual(double a, double b)
{
    double eps = 1e-10;
    if (fabs(a) <= eps || fabs(b) == eps)
        if (fabs(a - b) <= 2*eps)
            return true;
        else
            return false;
    else
        if (fabs(a - b) <= 2*eps*fabs(a) && fabs(a - b) <= 2*eps*fabs(b))
            return true;
        else
            return false;
    return false;
}

arma::mat faceSplitting(arma::mat A, const arma::mat &B)
{
    #pragma omp parallel for
    for (size_t i = 0; i < A.n_rows; i++)
        A.row(i) *= B(i);
    return A;
}

arma::mat operator^=(arma::mat A, const arma::mat &B)
{
    #pragma omp parallel for
    for (size_t i = 0; i < A.n_rows; i++)
        A.row(i) *= B(i);
    return A;
}

arma::mat operator&=(arma::mat A, const arma::mat &B)
{
    #pragma omp parallel for
    for (size_t i = 0; i < A.n_rows; i++)
        A.row(i) /= B(i);
    return A;
}

arma::mat operator%=(arma::mat A, const arma::vec &b)
{
    #pragma omp parallel for
    for (size_t i = 0; i < b.size(); i++)
        A.row(i) *= b(i);
    return A;
}

arma::mat operator/=(arma::mat A, const arma::vec &b)
{
    #pragma omp parallel for
    for (size_t i = 0; i < b.size(); i++)
        A.row(i) /= b(i);
    return A;
}