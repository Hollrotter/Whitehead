#pragma once
#include "BasisFunction.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

class ChebyshevU : public BasisFunction
{
public:
    ChebyshevU(size_t _n) : BasisFunction{_n}
    {
        xi =-cos(arma::datum::pi*(arma::regspace(0, n-1)+1)/(n+1));
    }
    virtual inline constexpr double constant() override
    {
        return 1;
    }
    virtual inline double linear(double x) override
    {
        return 2*x;
    }
    virtual inline constexpr double constantDerivative() override
    {
        return 0;
    }
    virtual inline constexpr double linearDerivative() override
    {
        return 2;
    }
    virtual inline void next(size_t k, double x, double &U, double &Up1) override
    {
        Up1 = boost::math::chebyshev_next(x, U, Up1);
    }
    virtual inline void nextDerivative(size_t k, double x, double &U, double &, double &dUp1) override
    {
        dUp1 = (k == 0) ? 8*x : 2*(k+2)*U + dUp1;
    }
    virtual inline double left(size_t k) override
    {
        return k%2==0 ? k+1 : -1.*(k+1);
    }
    virtual inline double right(size_t k) override
    {
        return k+1;
    }
    virtual inline double leftDerivative(size_t k) override
    {
        return k%2==0 ? -1.*k*(k+1)*(k+2)/3 : k*(k+1)*(k+2)/3;
    }
    virtual inline double rightDerivative(size_t k) override
    {
        return k*(k+1)*(k+2)/3;
    }
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) override
    {
        arma::vec c(k/2+1, arma::fill::none), f(k+1, arma::fill::zeros), df(k, arma::fill::zeros);
        for (size_t s = 0; s <= k/2; s++)
            c(s) = pow(-1, s)*pow(2, k-2*s)*bi(k-s, s);
        for (size_t s = 0; s <= k/2; s++)
        {
            arma::vec F(k+1, arma::fill::zeros);
            for (size_t t = 0; t <= k-2*s; t++)
                F(t) += bi(k-2*s, t) * pow(x, k-2*s-t);
            f += c(s) * F;
        }
        if (k > 0)
            for (size_t s = 0; s <= (k-1)/2; s++)
            {
                arma::vec F(k, arma::fill::zeros);
                for (size_t t = 0; t <= k-2*s-1; t++)
                    F(t) += bi(k-2*s-1, t)*pow(x, k-2*s-1-t);
                df += (k-2*s)*c(s) * F;
            }
        return std::make_pair(f, df);
    }
    friend class Wing;
};