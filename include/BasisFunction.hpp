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
    virtual inline constexpr double constant() const = 0;
    virtual inline double linear(double) const = 0;
    virtual inline constexpr double constantDerivative() const = 0;
    virtual inline constexpr double linearDerivative() const = 0;
    virtual inline void next(size_t, double, double&, double&) const = 0;
    virtual inline void nextDerivative(size_t, double, double&, double&, double&) const = 0;
    virtual inline double left(size_t) const = 0;
    virtual inline double right(size_t) const = 0;
    virtual inline double leftDerivative(size_t) const = 0;
    virtual inline double rightDerivative(size_t) const = 0;
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t, double, arma::umat&) const = 0;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t, size_t, arma::umat&, arma::mat&) const = 0;
    virtual inline double weightFunction(double) const = 0;
    virtual inline double weightFunctionDerivative(double) const = 0;
};