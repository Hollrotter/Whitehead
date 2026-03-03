#pragma once
#include "B_Spline.hpp"

class Splinefit
{
    size_t order = 3; // Order of the B-Splines
    arma::vec h; // Factors obtained from least-squares fitting
    size_t M; // Number of segments
    B_Spline B; // B-Spline
    arma::mat Q; // Orthogonal matrix
    arma::mat R; // Upper triangular matrix
    void setQR(const arma::vec);
public:
    Splinefit() = default;
    // Constructor for splinefitting (fitting itself is not done)
    Splinefit(arma::vec, size_t);
    // Constructor for splinefitting with immediate fitting
    Splinefit(arma::vec, arma::vec, size_t);
    // Fitting the B-Spline to given data
    void fit(const arma::vec y)
    {
        h = solve(trimatu(R), Q*y);
    }
    // Evaluate splinefitting at given value x returning a scalar z
    double operator()(const double);
    // Evaluate splinefitting at given vector x returning a vector z
    arma::vec operator()(const arma::vec);
    // Compute derivative of splinefitting at given value x returning a scalar dz/dx
    double diff(const double);
    // Compute derivaitve of splinefitting at given vector x returing a vector dz/dx
    arma::vec diff(const arma::vec);
};