#include "Chebyshev.hpp"

/**
 * @brief 
 * 
 * @param k 
 * @param x 
 * @return arma::vec 
 */
arma::vec Chebyshev::Polynomial(const size_t k, const arma::vec x)
{
	if (k == 0)
		return arma::ones(x.size());
	if (k == 1)
		return x;
	if (k <= 70)
	{
		arma::vec T_k;
		arma::vec T_k_2 = arma::ones(x.size());
		arma::vec T_k_1 = x;
		for (size_t j = 2; j <= k; j++)
		{
			T_k = 2*x%T_k_1 - T_k_2;
			T_k_2 = T_k_1;
			T_k_1 = T_k;
		}
		return T_k;
	}
	else
		return cos(k*acos(x));
}

/**
 * @brief 
 * 
 * @param u Vector of amplitudes (needed for the length in current implementation).
 * @param x Vector of values where the Polynomial will be evaluated.
 * @return arma::mat 
 */
arma::mat Chebyshev::Polynomial(const arma::vec u, const arma::vec x)
{
	arma::mat T(x.size(), u.size());
    #pragma omp parallel for
	for (size_t k = 0; k < u.size(); k++)
		T.col(k) = Polynomial(k, x);
	return T;
}