#pragma once
#include "BasisFunction.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

class ChebyshevU : public BasisFunction
{
public:
    ChebyshevU(size_t _n, size_t _m);
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
        return 2*x;
    }
    virtual inline arma::vec linear(arma::vec x) const override
    {
        return 2*x;
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
        return 2;
    }
    virtual inline arma::vec linearDerivative(size_t k) const override
    {
        return 2*arma::ones(k);
    }
    virtual inline void next(size_t _, double x, double &U, double &Up1) const override
    {
        Up1 = boost::math::chebyshev_next(x, U, Up1);
    }
    virtual inline void next(size_t _, arma::vec x, arma::vec &U, arma::vec &Up1) const override
    {
        Up1 = 2*x%U - Up1;
    }
    virtual inline void nextDerivative(size_t k, double x, double &U, double &_, double &dUp1) const override
    {
        dUp1 = (k == 0) ? 8*x : 2*(k+2)*U + dUp1;
    }
    virtual inline void nextDerivative(size_t k, arma::vec x, arma::vec &U, arma::vec &_, arma::vec &dUp1) const override
    {
        dUp1 = (k == 0) ? arma::vec(8*x) : arma::vec(2*(k+2)*U + dUp1);
    }
    virtual inline double left(size_t k) const override
    {
        return k%2==0 ? k+1 : -1.*(k+1);
    }
    virtual inline double right(size_t k) const override
    {
        return k+1;
    }
    virtual inline double leftDerivative(size_t k) const override
    {
        return k%2==0 ? -1.*k*(k+1)*(k+2)/3 : k*(k+1)*(k+2)/3;
    }
    virtual inline double rightDerivative(size_t k) const override
    {
        return k*(k+1)*(k+2)/3;
    }
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x) const override;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t i, arma::mat &bi_05) const override;
    virtual inline double weightFunction(double x) const override
    {
        return sqrt(1 - x*x);
    }
    virtual inline arma::vec weightFunction(arma::vec x) const override
    {
        return sqrt(1 - x%x);
    }
    virtual inline double weightFunctionDerivative(double x) const override
    {
        return -x/sqrt(1 - x*x);
    }
    virtual inline arma::vec weightFunctionDerivative(arma::vec x) const override
    {
        return -x/sqrt(1 - x%x);
    }
    friend class Wing;
};