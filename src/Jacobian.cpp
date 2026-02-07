#include "Metric.hpp"

/**
 * @brief 
 * 
 * @param y1 Matrix of values of first variable.
 * @param y2 Matrix of values of second variable.
 * @param D1 Derivative Matrix of differentiation in 1-direction.
 * @param D2 Derivative Matrix of differentiation in 2-direction.
 * @return arma::field<arma::mat> 
 */
arma::field<arma::mat> jacobian(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2)
{
	arma::field<arma::mat> J(2, 2);

	try
	{
		J(0, 0) = D1 * y1;
	}
	catch(const std::exception& e)
	{
		std::println("D1({} x {}) and y1({} x {}) mismatch for jacobian!", D1.n_rows, D1.n_cols, y1.n_rows, y1.n_cols);
		exit(EXIT_FAILURE);
	}
	try
	{
		J(0, 1) = y1 * D2.t();
	}
	catch(const std::exception& e)
	{
		std::println("y1({} x {}) and D2.t()({} x {}) mismatch for jacobian!", y1.n_rows, y1.n_cols, D2.n_cols, D2.n_rows);
		exit(EXIT_FAILURE);
	}
	try
	{
		J(1, 0) = D1 * y2;
	}
	catch(const std::exception& e)
	{
		std::println("D1({} x {}) and y2({} x {}) mismatch for jacobian!", D1.n_rows, D1.n_cols, y2.n_rows, y2.n_cols);
		exit(EXIT_FAILURE);
	}
	try
	{
		J(1, 1) = y2 * D2.t();
	}
	catch(const std::exception& e)
	{
		std::println("y2({} x {}) and D2.t()({} x {}) mismatch for jacobian!", y2.n_rows, y2.n_cols, D2.n_cols, D2.n_rows);
		exit(EXIT_FAILURE);
	}
	return J;
}