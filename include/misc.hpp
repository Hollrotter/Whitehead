#pragma once
#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <armadillo>
#include <print>
#include <omp.h>

template <typename T>
concept Number = std::integral<T> || std::floating_point<T>;
arma::mat faceSplitting(arma::mat A, const arma::mat &B);
inline arma::mat QLM(const arma::mat &A, arma::mat B)
{
    return diagmat(vectorise(A))*B;
}
arma::mat operator^=(arma::mat A, const arma::mat &B);
arma::mat operator&=(arma::mat A, const arma::mat &B);
arma::mat operator%=(arma::mat A, const arma::vec &b);
arma::mat operator/=(arma::mat A, const arma::vec &b);
inline arma::mat operator^(const arma::mat &A, arma::mat B)
{
    return B ^= A;
}
inline arma::mat operator&(arma::mat A, const arma::mat &B)
{
    return A &= B;
}
inline arma::mat operator%(arma::mat A, const arma::vec &b)
{
    return A %= b;
}
inline arma::mat operator%(const arma::vec &b, arma::mat A)
{
    return A %= b;
}
inline arma::mat operator/(arma::mat A, const arma::vec &b)
{
    return A /= b;
}