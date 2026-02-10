#include "Chebyshev.hpp"

/**
 * @brief 
 * 
 * @param x Vector of nodes forming the Gauss-Lobatto Nodes.
 * @param u Function values at the given nodes x.
 * @return arma::vec 
 */
arma::vec Chebyshev::DiscreteChebyshevTransform(const arma::vec x, const arma::vec u)
{
	size_t N = x.size();
    if (u.size() != N)
    {
        std::println("u must have the same size as x for DiscreteChebyshevTransform!");
        exit(EXIT_FAILURE);
    }
	arma::vec u_hat(N, arma::fill::zeros);
	const arma::mat T = Polynomial(u_hat,-x);
    #pragma omp parallel for
	for (size_t k = 0; k < N; k++)
	{
		for (size_t j = 0; j < N; j++)
			u_hat(k) += u(j)*T(k, j)*w(j, N-1);
		u_hat(k) /= gamma(k, N-1);
	}
	return u_hat;
}

/**
 * @brief 
 * 
 * @param x Vector of nodes forming the Gauss-Lobatto Nodes.
 * @param y Vector of nodes forming the Gauss-Lobatto Nodes.
 * @param u Function values at the given grid defined by x and y
 * @return arma::mat 
 */
arma::mat Chebyshev::DiscreteChebyshevTransform(const arma::vec x, const arma::vec y, const arma::mat u)
{
	size_t N = x.size();
	size_t M = y.size();
    if (u.n_rows != N || u.n_cols != M)
    {
        std::println("u must have the shape of x and y for DiscreteChebyshevTransform!");
        exit(EXIT_FAILURE);
    }
	arma::mat u_hat(N, M, arma::fill::zeros);
	const arma::mat Tx = Polynomial(u_hat.col(0),-x);
	const arma::mat Ty = Polynomial(u_hat.row(0).t(),-y);
    #pragma omp parallel for
	for (size_t n = 0; n < N; n++)
		for (size_t m = 0; m < M; m++)
		{
			for (size_t i = 0; i < N; i++)
				for (size_t j = 0; j < M; j++)
					u_hat(n, m) += u(i, j)*Tx(n, i)*Ty(m, j)*w(i, N-1)*w(j, M-1);
			u_hat(n, m) /= gamma(n, N-1)*gamma(m, M-1);
		}
	return u_hat;
}

/**
 * @brief 
 * 
 * @param u_hat Amplitudes calculated by DiscreteChebyshevTransform(x, u)
 * @return arma::vec 
 */
arma::vec Chebyshev::ChebyshevDerivativeCoefficients(const arma::vec u_hat)
{
	size_t n = u_hat.size();
	arma::vec u_bar(n);
	u_bar(n-1) = 0;
	u_bar(n-2) = 2*n*u_hat(n-1);
    #pragma omp parallel for
	for (size_t k = n-3; k > 0; k--)
		u_bar(k) = 2*(k+1)*u_hat(k+1) + u_bar(k+2);
	u_bar(0) = u_hat(1) + u_bar(2)/2;
	return u_bar;
}

/**
 * @brief 
 * 
 * @param u_hat 
 * @return arma::cube 
 */
arma::cube Chebyshev::ChebyshevDerivativeCoefficients(const arma::mat u_hat)
{
	size_t n = u_hat.n_rows;
	size_t m = u_hat.n_cols;
	arma::cube u_bar(n, m, 2);
	u_bar.slice(0).row(n-1).zeros();
	u_bar.slice(0).row(n-2) = 2*n*u_hat.row(n-1);
	u_bar.slice(1).col(m-1).zeros();
	u_bar.slice(1).col(m-2) = 2*m*u_hat.col(m-1);
    #pragma omp parallel for
	for (size_t k = n-3; k > 0; k--)
		u_bar.slice(0).row(k) = 2*(k+1)*u_hat.row(k+1) + u_bar.slice(0).row(k+2);
    #pragma omp parallel for
	for (size_t k = m-3; k > 0; k--)
		u_bar.slice(1).col(k) = 2*(k+1)*u_hat.col(k+1) + u_bar.slice(1).col(k+2);
	u_bar.slice(0).row(0) = u_hat.row(1) + u_bar.slice(0).row(2)/2;
	u_bar.slice(1).col(0) = u_hat.col(1) + u_bar.slice(1).col(2)/2;
	return u_bar;
}

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