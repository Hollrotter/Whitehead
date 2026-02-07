#include "Metric.hpp"

/**
 * @brief 
 * 
 * @param y1 
 * @param y2 
 * @param D1 
 * @param D2 
 * @return std::tuple<arma::cube, arma::cube> 
 */
std::tuple<arma::cube, arma::cube> Metric(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2)
{
	arma::cube g_c = MetricCo(y1, y2, D1, D2);
	arma::cube gc  = MetricContra(g_c);
	return {g_c, gc};
}

/**
 * @brief 
 * 
 * @param J 
 * @return std::tuple<arma::cube, arma::cube> 
 */
std::tuple<arma::cube, arma::cube> Metric(const arma::field<arma::mat> J)
{
	arma::cube g_c = MetricCo(J);
	arma::cube gc  = MetricContra(g_c);
	return {g_c, gc};
}

/**
 * @brief 
 * 
 * @param g 
 * @return arma::cube 
 */
arma::cube MetricContra(arma::cube &g)
{
	arma::cube gc(g.n_rows, g.n_cols, 3);

	#pragma omp parallel for
	for (size_t i = 0; i < g.n_rows; i++)
		for (size_t j = 0; j < g.n_cols; j++)
			gc.tube(i, j) = arma::vec{g(i, j, 2),-g(i, j, 1), g(i, j, 0)}/(g(i, j, 0)*g(i, j, 2)-pow(g(i, j, 1), 2));
	return gc;
}