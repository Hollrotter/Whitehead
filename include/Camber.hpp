#pragma once
#include "enums.hpp"
#include "Splinefit.hpp"

class Camber
{
    CamberType camberType; //Camber type (none, function or B_Spline)
    std::function<double(double)>  f; // Analytical function for the camber
    std::function<double(double)> df; // Analytical function for the camber slope
    Splinefit s; // Splinefitting for the initial contour of the camber
public:
    Camber() : camberType(CamberType::none) {}
    explicit Camber(std::function<double(double)> dF) : camberType(CamberType::function), df(dF) {}
    Camber(std::function<double(double)> F, std::function<double(double)> dF) : camberType(CamberType::function), f(F), df(dF) {}
    explicit Camber(const Splinefit &S) : camberType(CamberType::b_spline), s(S) {}
    // Compute z coordinate at given x value
    double operator() (double);
    // Compute z coordinates at given vector x
    arma::vec operator() (arma::vec);
    // Compute derivative at given x value
    double diff(double);
    // Compute derivatives at given vector x
    arma::vec diff(arma::vec);
};