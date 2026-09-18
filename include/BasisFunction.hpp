#pragma once
#include "misc.hpp"

class BasisFunction
{
protected:
    size_t n;
    arma::vec xi;
public:
    BasisFunction() = default;
    BasisFunction(size_t _n) : n(_n) {};
    virtual inline constexpr double constant() = 0;
    virtual inline double linear(double) = 0;
    virtual inline constexpr double constantDerivative() = 0;
    virtual inline constexpr double linearDerivative() = 0;
    virtual inline void next(size_t, double, double&, double&) = 0;
    virtual inline void nextDerivative(size_t, double, double&, double&, double&) = 0;
    virtual inline double left(size_t) = 0;
    virtual inline double right(size_t) = 0;
    virtual inline double leftDerivative(size_t) = 0;
    virtual inline double rightDerivative(size_t) = 0;
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t, double, arma::umat&) = 0;
};