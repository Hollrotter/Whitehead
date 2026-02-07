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
arma::cube Christoffel(arma::cube &g_c, arma::cube &gc, arma::mat &D1, arma::mat &D2)
{
	arma::mat dg11d1, dg11d2, dg12d1, dg12d2, dg22d1, dg22d2;
	try
	{
		dg11d1 = D1*g_c.slice(0);
	}
	catch(const std::exception& e)
	{
		std::println("D1({} x {}) and g_11({} x {}) mismatch for christoffel!", D1.n_rows, D1.n_cols, g_c.slice(0).n_rows, g_c.slice(0).n_cols);
	}
	try
	{
		dg11d2 = g_c.slice(0)*D2.t();
	}
	catch(const std::exception& e)
	{
		std::println("g_11({} x {}) and D2.t()({} x {}) mismatch for christoffel!", g_c.slice(0).n_rows, g_c.slice(0).n_cols, D2.n_cols, D2.n_rows);
	}
	try
	{
		dg12d1 = D1*g_c.slice(1);
	}
	catch(const std::exception& e)
	{
		std::println("D1({} x {}) and g_12({} x {}) mismatch for christoffel!", D1.n_rows, D1.n_cols, g_c.slice(1).n_rows, g_c.slice(1).n_cols);
	}
	try
	{
		dg12d2 = g_c.slice(1)*D2.t();
	}
	catch(const std::exception& e)
	{
		std::println("g_12({} x {}) and D2.t()({} x {}) mismatch for christoffel!", g_c.slice(1).n_rows, g_c.slice(1).n_cols, D2.n_cols, D2.n_rows);
	}
	try
	{
		dg22d1 = D1*g_c.slice(2);
	}
	catch(const std::exception& e)
	{
		std::println("D1({} x {}) and g_22({} x {}) mismatch for christoffel!", D1.n_rows, D1.n_cols, g_c.slice(2).n_rows, g_c.slice(2).n_cols);
	}
	try
	{
		dg22d2 = g_c.slice(2)*D2.t();
	}
	catch(const std::exception& e)
	{
		std::println("g_22({} x {}) and D2.t()({} x {}) mismatch for christoffel!", g_c.slice(2).n_rows, g_c.slice(2).n_cols, D2.n_cols, D2.n_rows);
	}
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