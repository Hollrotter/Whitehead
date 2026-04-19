#pragma once
#include "Lagrange.hpp"

arma::field<arma::mat> Jacobian(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2);
arma::field<arma::mat> Jacobian(std::array<Lagrange::CurveInterpolant*, 4> chi);
arma::field<arma::mat> Jacobian(arma::vec &x1, arma::vec &x2, std::array<Lagrange::CurveInterpolant*, 4> chi);
std::tuple<arma::cube, arma::cube> Metric(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2);
std::tuple<arma::cube, arma::cube> Metric(const std::array<Lagrange::CurveInterpolant*, 4> chi);
std::tuple<arma::cube, arma::cube> Metric(const arma::field<arma::mat> J);
arma::cube MetricCo(arma::mat &y1, arma::mat &y2, arma::mat &D1, arma::mat &D2);
arma::cube MetricCo(const std::array<Lagrange::CurveInterpolant*, 4> chi);
arma::cube MetricCo(const arma::field<arma::mat> J);
arma::cube MetricContra(const arma::cube &g);
arma::cube Christoffel(const arma::cube &g_c, const arma::cube &gc, const arma::mat &D1, const arma::mat &D2);