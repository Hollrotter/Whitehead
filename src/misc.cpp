#include "misc.hpp"

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