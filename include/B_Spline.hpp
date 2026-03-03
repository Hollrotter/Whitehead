#pragma once
#include "misc.hpp"

class B_Spline
{
    size_t order; // Polynomial order of the Spline
    arma::vec X; // Nodes of the Spline segments
public:
    B_Spline() = default;
    // Constructor setting the segment nodes and the order of the spline
    B_Spline(arma::vec, size_t);
    arma::vec operator()(const double);
    // Evaluate spline basis i at given value x with order of spline
    double operator()(const double x, const size_t i)
    {
        return operator()(x, i, order);
    }
    // Evaluate spline basis i at given value x with arbitrary order k
    double operator()(const double, const size_t, const size_t);
    arma::mat operator()(const arma::vec);
    // Evaluate spline basis i at given vector x with order of spline
    arma::vec operator()(const arma::vec x, const size_t i)
    {
        return operator()(x, i, order);
    }
    // Evaluate spline basis i at given vector x with arbitrary order k
    arma::vec operator()(const arma::vec, const size_t, const size_t);
    arma::vec diff(const double);
    // Compute derivative from spline approximation
    double diff(const double x, const size_t i)
    {
        return diff(x, i, order);
    }
    // Compute derivative from spline approximation
    double diff(const double, const size_t, const size_t);
    // Compute derivative from spline approximation
    arma::vec diff(const arma::vec x, const size_t i)
    {
        return diff(x, i, order);
    }
    // Compute derivative from spline approximation
    arma::vec diff(const arma::vec, const size_t, const size_t);
    double integrate(const size_t i, const size_t j);
    double integrate(const size_t i, const size_t j, const size_t k);
    double integrate(const double x1, const double x2, const size_t j);
    double integrate(const double x1, const double x2, const size_t j, const size_t k);
    double integrate_x(const double x1, const double x2, const size_t j);
};