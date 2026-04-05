#pragma once
#include "Lagrange.hpp"

arma::mat ChebyshevInterpolation(arma::vec x1, arma::vec x2, arma::mat xs, arma::mat ys, arma::mat xss, arma::mat yss);
arma::mat LagrangeInterpolation(arma::vec x1, arma::vec x2, arma::mat xs, arma::mat ys, arma::mat xss, arma::mat yss);