#pragma once
#include "BasisFunction.hpp"
#include <boost/math/special_functions/chebyshev.hpp>

class ChebyshevU : public BasisFunction
{
public:
    ChebyshevU(size_t _n, size_t _m) : BasisFunction{_n, _m}
    {
        xi =-cos(arma::datum::pi*(arma::regspace(0, n-1)+1)/(n+1));
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
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) const override
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
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t i, size_t k, arma::umat &bi, arma::mat &bi_05) const override
    {
        arma::vec c(k+1), d(k+1);
        for (size_t ii = 0; ii <= k; ii++)
            for (size_t jj = 0; jj <= ii/2; jj++)
                c(ii) += pow(-1, ii-jj)*bi_05(1, ii-jj)*bi(ii-jj, jj)*pow(2*xi(i), ii-2*jj)/pow(1-pow(xi(i), 2), ii-jj);
        c *= sqrt(1-pow(xi(i), 2));
        for (size_t ii = 0; ii < k; ii++)
            for (size_t jj = 0; jj <= ii/2; jj++)
            {
                double con = pow(-1, ii-jj)*bi(ii-jj, jj)*pow(2*xi(i), ii-2*jj)/pow(1-pow(xi(i), 2), ii-jj);
                d(ii)   += con * bi_05(0, ii-jj) * xi(i);
                d(ii+1) += con * bi_05(0, ii-jj);
            }
        d /=-sqrt(1-pow(xi(i), 2));
        return std::make_pair(c, d);
    }
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