#pragma once
#include "BasisFunction.hpp"
#include "Chebyshev.hpp"
#include "fastgl.h"

class ChebyshevT : public BasisFunction
{
public:
    ChebyshevT(size_t _n, size_t _);
    virtual inline constexpr double constant() const override
    {
        return 1;
    }
    virtual inline arma::vec constant(size_t k) const override
    {
        return arma::ones(k);
    }
    virtual inline double linear(double x) const override
    {
        return x;
    }
    virtual inline arma::vec linear(arma::vec x) const override
    {
        return x;
    }
    virtual inline constexpr double constantDerivative() const override
    {
        return 0;
    }
    virtual inline arma::vec constantDerivative(size_t k) const override
    {
        return arma::zeros(k);
    }
    virtual inline constexpr double linearDerivative() const override
    {
        return 1;
    }
    virtual inline arma::vec linearDerivative(size_t k) const override
    {
        return arma::ones(k);
    }
    virtual inline void next(size_t _, double x, double &T, double &Tp1) const override
    {
        Tp1 = boost::math::chebyshev_next(x, T, Tp1);
    }
    virtual inline void next(size_t k, arma::vec x, arma::mat &T) const override
    {
        T.col(k+1) = 2*x%T.col(k) - T.col(k-1);
    }
    virtual inline void nextDerivative(size_t k, double x, double &T, double &_, double &dTp1) const override
    {
        dTp1 = (k == 0) ? 4*x : (k+2)*(2*T + dTp1/k);
    }
    virtual inline void nextDerivative(size_t k, arma::vec _, arma::mat &T, arma::mat &dT) const override
    {
        dT.col(k+1) = arma::vec((k+2)*(2*T.col(k) + dT.col(k-1)/k));
    }
    virtual inline double left(size_t k) const override
    {
        return k%2==0 ? 1 : -1;
    }
    virtual inline double right(size_t _) const override
    {
        return 1;
    }
    virtual inline double leftDerivative(size_t k) const override
    {
        return k%2==0 ? -1.*k*k : k*k;
    }
    virtual inline double rightDerivative(size_t k) const override
    {
        return k*k;
    }
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x) const override;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t _, arma::mat &bi_05) const override;
    virtual inline double weightFunction(double _) const override
    {
        return 1;
    }
    virtual inline arma::vec weightFunction(arma::vec x) const override
    {
        return arma::ones(x.size());
    }
    virtual inline double weightFunctionDerivative(double _) const override
    {
        return 0;
    }
    virtual inline arma::vec weightFunctionDerivative(arma::vec x) const override
    {
        return arma::zeros(x.size());
    }
    friend class Wing;
};