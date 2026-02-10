#pragma once
#include "misc.hpp"
#include "enums.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

namespace Chebyshev
{
    // Evaluates the Chebyshev Polynomial of order k at values x.
    arma::vec Polynomial(const size_t k, const arma::vec x);
    // Evaluates all the Chebyshev Polynomials up to order len(u) at values x.
    arma::mat Polynomial(const arma::vec u, const arma::vec x);
    // Calculates the amplitudes for Discrete-Chebyshev-Transform in 1D
    arma::vec DiscreteChebyshevTransform(const arma::vec x, const arma::vec u);
    // Calculates the amplitudes for Discrete-Chebyshev-Transform in 2D
    arma::mat DiscreteChebyshevTransform(const arma::vec x, const arma::vec y, const arma::mat u);
    arma::vec  ChebyshevDerivativeCoefficients(const arma::vec u_hat);
    arma::cube ChebyshevDerivativeCoefficients(const arma::mat u_hat);
    arma::cx_vec FFFTEO(const arma::vec f);
    std::pair<arma::vec, arma::vec> forwardRealFFT(const arma::vec x);
    arma::vec fastCosineTransform(const arma::vec f, const std::string s);
    arma::vec fastChebyshevTransform(const arma::vec f, const std::string s);
    arma::mat fastChebyshevTransform(const arma::mat f, const std::string s);
    inline double c(const size_t k)
    {
        return (k == 0) ? 2 : 1;
    }
    inline double gamma(const size_t k, const size_t N)
    {
        return (k == N) ? arma::datum::pi : arma::datum::pi/2*c(k);
    }
    inline double w(const size_t j, const size_t N)
    {
        return (j == 0 || j == N) ? arma::datum::pi/2/N : arma::datum::pi/N;
    }
    // Generates the Gauss-Lobatto-Nodes with length n
    inline arma::vec gaussLobatto(const size_t n)
    {
        return -cos(arma::datum::pi*arma::regspace(0, n-1)/(n-1));
    }
    // Calculates the derivative Matrix for Gauss-Lobatto-Nodes x for given derivative.
    arma::mat derivativeMatrix(const arma::vec x, const Derivative derivative);
}