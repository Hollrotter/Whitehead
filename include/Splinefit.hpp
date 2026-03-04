#pragma once
#include "B_Spline.hpp"

class Splinefit
{
    size_t order = 3; // Order of the B-Splines
    arma::vec h; // Factors obtained from least-squares fitting
    size_t M; // Number of segments in 1-direction
    size_t N; // Number of segments in 2-direction
    B_Spline B1; // B-Spline in 1-direction
    B_Spline B2; // B-Spline in 2-direction
    arma::mat Q; // Orthogonal matrix
    arma::mat R; // Upper triangular matrix
    void setQR(const arma::vec);
    void setQR(const arma::vec, const arma::vec);
public:
    Splinefit() = default;
    // Constructor for 1D splinefitting (fitting itself is not done)
    Splinefit(arma::vec, size_t);
    // Constructor for 1D splinefitting with immediate fitting
    Splinefit(arma::vec, arma::vec, size_t);
    // Constructor for 2D splinefitting (fitting itself is not done)
    Splinefit(arma::vec, arma::vec, size_t, size_t);
    // Constructor for 2D splinefitting with immediate fitting
    Splinefit(arma::vec, arma::vec, arma::mat, size_t, size_t);
    // Fitting the B-Spline to given 1D-data
    void fit(const arma::vec y)
    {
        h = solve(trimatu(R), Q*y);
    }
    // Fitting the B-Spline to given 2D-data
    void fit(const arma::mat z)
    {
        h = solve(trimatu(R), Q*vectorise(z));
    }
    // Evaluate 1D splinefitting at given value x returning a scalar z
    double operator()(const double);
    // Evaluate 1D splinefitting at given vector x returning a vector z
    arma::vec operator()(const arma::vec);
    // Evaluate 2D splinefitting at given values x and y returning a scalar z
    double operator()(const double, const double);
    // Evaluate 2D splinefitting at given vectors x and y returning a matrix z
    arma::mat operator()(const arma::vec, const arma::vec);
    // Compute derivative of 1D splinefitting at given value x returning a scalar dz/dx
    double diff(const double);
    // Compute derivaitve of 1D splinefitting at given vector x returing a vector dz/dx
    arma::vec diff(const arma::vec);
    // Compute derivative of 2D splinefitting at given values x and y returning gradient vector
    arma::vec diff(const double, const double);
    // Compute derivative of 2D splinefitting at given vectors x and y returning gradient cube
    arma::cube diff(const arma::vec, const arma::vec);
};