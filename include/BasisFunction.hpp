#pragma once
#include "gaujac.hpp"

class BasisFunction
{
protected:
    size_t n;
    size_t m;
    arma::vec xi;
    arma::vec xg;
    arma::vec wg;
public:
    BasisFunction() = default;
    BasisFunction(size_t _n, size_t _m) : n(_n), m(_m) {};
    virtual inline constexpr double constant() const = 0;
    virtual inline arma::vec constant(size_t) const = 0;
    virtual inline double linear(double) const = 0;
    virtual inline arma::vec linear(arma::vec) const = 0;
    virtual inline constexpr double constantDerivative() const = 0;
    virtual inline arma::vec constantDerivative(size_t) const = 0;
    virtual inline constexpr double linearDerivative() const = 0;
    virtual inline arma::vec linearDerivative(size_t) const = 0;
    virtual inline void next(size_t, double, double&, double&) const = 0;
    virtual inline void next(size_t, arma::vec, arma::mat&) const = 0;
    virtual inline void nextDerivative(size_t, double, double&, double&, double&) const = 0;
    virtual inline void nextDerivative(size_t, arma::vec, arma::mat&, arma::mat&) const = 0;
    virtual inline double left(size_t) const = 0;
    virtual inline double right(size_t) const = 0;
    virtual inline double leftDerivative(size_t) const = 0;
    virtual inline double rightDerivative(size_t) const = 0;
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t, double) const = 0;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t, arma::mat&) const = 0;
    virtual inline double weightFunction(double) const = 0;
    virtual inline arma::vec weightFunction(arma::vec) const = 0;
    virtual inline double weightFunctionDerivative(double) const = 0;
    virtual inline arma::vec weightFunctionDerivative(arma::vec) const = 0;
    friend class Wing;
};