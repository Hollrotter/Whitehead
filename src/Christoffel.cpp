#include "Metric.hpp"

/**
 * @brief 
 * 
 * @param g_c 
 * @param gc 
 * @param D1 
 * @param D2 
 * @return arma::cube 
 */
arma::cube Christoffel(const arma::cube &g_c, const arma::cube &gc, const arma::mat &D1, const arma::mat &D2)
{
	arma::cube dgd1, dgd2;
	try
	{
		dgd1 = arma::cubemul(D1, g_c);
	}
	catch(const std::exception& e)
	{
		std::println("Size mismatch in calculation of Christoffel Symbols!");
	}
	try
	{
		dgd2 = arma::cubemul(g_c, D2.t());
	}
	catch(const std::exception& e)
	{
		std::println("Size mismatch in calculation of Christoffel Symbols!");
	}
	arma::mat dg11d1 = dgd1.slice(0);
	arma::mat dg12d1 = dgd1.slice(1);
	arma::mat dg22d1 = dgd1.slice(2);
	arma::mat dg11d2 = dgd2.slice(0);
	arma::mat dg12d2 = dgd2.slice(1);
	arma::mat dg22d2 = dgd2.slice(2);
	arma::cube gam(g_c.n_rows, g_c.n_cols, 6);
	arma::mat g11 = gc.slice(0);
	arma::mat g12 = gc.slice(1);
	arma::mat g22 = gc.slice(2);
	gam.slice(0) = (g11%dg11d1 + g12%(2*dg12d1 - dg11d2))/2;
	gam.slice(1) = (g11%dg11d2 + g12%dg22d1)/2;
	gam.slice(2) = (g11%(2*dg12d2 - dg22d1) + g12%dg22d2)/2;
	gam.slice(3) = (g12%dg11d1 + g22%(2*dg12d1 - dg11d2))/2;
	gam.slice(4) = (g12%dg11d2 + g22%dg22d1)/2;
	gam.slice(5) = (g12%(2*dg12d2 - dg22d1) + g22%dg22d2)/2;
	return gam;
}