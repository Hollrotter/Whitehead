#include "Chebyshev.hpp"

/**
 * @brief 
 * 
 * @param x Vector of nodes forming the Gauss-Lobatto Nodes.
 * @param derivative Derivative to be calculated (First or Second are implemented).
 * @return arma::mat 
 */
arma::mat Chebyshev::derivativeMatrix(const arma::vec x, const Derivative derivative)
{
    size_t n = x.size()-1;
    arma::mat D(n+1, n+1);
    auto cBar = [] (size_t i, size_t n) {return ((i == 0) || (i == n) ? 2. : 1.);};
    switch (derivative)
    {
        case Derivative::first:
            #pragma omp parallel for
            for (size_t i = 0; i <= n; i++)
                for (size_t j = 0; j <= n; j++)
                    if (i == 0 && j == 0)
                        D(i, j) =-(2*pow(n, 2) + 1)/6;
                    else if (i == n && j == n)
                        D(i, j) = (2*pow(n, 2) + 1)/6;
                    else if (i == j)
                        D(i, j) =-x(i)/(2*(1 - x(i)*x(i)));
                    else
                        D(i, j) = cBar(i, n)/cBar(j, n)*pow(-1, i+j)/(x(i) - x(j));
            return D;
        case Derivative::second:
            #pragma omp parallel for
            for (size_t i = 0; i <= n; i++)
                for (size_t j = 0; j <= n; j++)
                    if (i == 0 && j == 0 || i == n && j == n)
                        D(i, i) = (pow(n, 4) - 1)/15;
                    else if (i == j)
                        D(i, i) =-((n*n - 1)*(1 - x(i)*x(i)) + 3)/(3*pow(1 - x(i)*x(i), 2));
                    else if (i == 0)
                        D(0, j) = 2*pow(-1,   j)/3/cBar(j, n)*((2*n*n + 1)*(1 + x(j)) - 6)/pow(1 + x(j), 2);
                    else if (i == n)
                        D(n, j) = 2*pow(-1, j+n)/3/cBar(j, n)*((2*n*n + 1)*(1 - x(j)) - 6)/pow(1 - x(j), 2);
                    else
                        D(i, j) = pow(-1, i+j)/cBar(j, n)*(x(i)*x(i) + x(i)*x(j) - 2)/(1 - x(i)*x(i))/pow(x(i) - x(j), 2);
            return D;
    }
    std::unreachable();
}