#pragma once
#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <armadillo>
#include <omp.h>
#include "enums.hpp"

namespace Chebyshev
{
    // Generates the Gauss-Lobatto-Nodes with length n
    inline arma::vec gaussLobatto(const size_t n)
    {
        return -cos(arma::datum::pi*arma::regspace(0, n-1)/(n-1));
    }
    // Calculates the derivative Matrix for Gauss-Lobatto-Nodes x for given derivative.
    arma::mat derivativeMatrix(const arma::vec x, const Derivative derivative);
}