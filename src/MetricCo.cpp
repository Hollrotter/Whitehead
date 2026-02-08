#include "Metric.hpp"

/**
 * @brief 
 * 
 * @param y1 Matrix of values of first variable.
 * @param y2 Matrix of values of second variable.
 * @param D1 Derivative Matrix of differentiation in 1-direction.
 * @param D2 Derivative Matrix of differentiation in 2-direction.
 * @return arma::cube 
 */
arma::cube MetricCo(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2)
{
	arma::mat dxdx1 = D1 * y1;
	arma::mat dxdx2 = y1 * D2.t();
	arma::mat dydx1 = D1 * y2;
	arma::mat dydx2 = y2 * D2.t();

	arma::cube g_c(D1.n_rows, D2.n_rows, 3);
	#pragma omp parallel for
	for (size_t i = 0; i < D1.n_rows; i++)
		for (size_t j = 0; j < D2.n_rows; j++)
		{
			arma::vec g1 = {dxdx1(i, j), dydx1(i, j)};
			arma::vec g2 = {dxdx2(i, j), dydx2(i, j)};
			g_c.tube(i, j) = arma::vec{dot(g1, g1), dot(g1, g2), dot(g2, g2)};
		}
	return g_c;
}

arma::cube MetricCo(const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
	size_t nx = chi[1]->getNodes().size();
	size_t ny = chi[0]->getNodes().size();
	arma::cube g_c(nx, ny, 3);
	auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(chi);
	#pragma omp parallel for
	for (size_t i = 0; i < nx; i++)
		for (size_t j = 0; j < ny; j++)
		{
			arma::vec g1 = {dxdx1(i, j), dydx1(i, j)};
			arma::vec g2 = {dxdx2(i, j), dydx2(i, j)};
			g_c.tube(i, j) = arma::vec{dot(g1, g1), dot(g1, g2), dot(g2, g2)};
		}
	return g_c;
}

/**
 * @brief 
 * 
 * @param J Jacobian
 * @return arma::cube 
 */
arma::cube MetricCo(const arma::field<arma::mat> J)
{
	arma::cube g_c(J(0, 0).n_rows, J(0, 0).n_cols, 3);

	g_c.slice(0) = pow(J(0, 0), 2) + pow(J(1, 0), 2);
	g_c.slice(1) = J(0, 0)%J(0, 1) + J(1, 0)%J(1, 1);
	g_c.slice(2) = pow(J(0, 1), 2) + pow(J(1, 1), 2);

	return g_c;
}