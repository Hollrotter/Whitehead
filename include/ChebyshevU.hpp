#pragma once
#include "BasisFunction.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

class ChebyshevU : public BasisFunction
{
public:
    ChebyshevU(size_t _n, size_t _m) : BasisFunction{_n, _m}
    {
        xi =-cos(arma::datum::pi*(arma::regspace(0, n-1)+1)/(n+1));
        xg = xi;
        wg = arma::datum::pi/(n+1)*(1-xi%xi);
    }
    virtual inline constexpr double constant() const override
    {
        return 1;
    }
    virtual inline double linear(double x) const override
    {
        return 2*x;
    }
    virtual inline constexpr double constantDerivative() const override
    {
        return 0;
    }
    virtual inline constexpr double linearDerivative() const override
    {
        return 2;
    }
    virtual inline void next(size_t k, double x, double &U, double &Up1) const override
    {
        Up1 = boost::math::chebyshev_next(x, U, Up1);
    }
    virtual inline void nextDerivative(size_t k, double x, double &U, double &, double &dUp1) const override
    {
        dUp1 = (k == 0) ? 8*x : 2*(k+2)*U + dUp1;
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
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) const override;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t i, size_t k, arma::umat &bi, arma::mat &bi_05) const override;
    virtual inline double weightFunction(double x) const override
    {
        return sqrt(1 - x*x);
    }
    virtual inline double weightFunctionDerivative(double x) const override
    {
        return -x/sqrt(1 - x*x);
    }
    friend class Wing;
};