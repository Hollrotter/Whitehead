#pragma once
#include "BasisFunction.hpp"
#include "Chebyshev.hpp"

class ChebyshevT : public BasisFunction
{
public:
    ChebyshevT(size_t _n) : BasisFunction{_n}
    {
        xi = Chebyshev::gauss(n);
    }
    virtual inline constexpr double constant() const override
    {
        return 1;
    }
    virtual inline double linear(double x) const override
    {
        return x;
    }
    virtual inline constexpr double constantDerivative() const override
    {
        return 0;
    }
    virtual inline constexpr double linearDerivative() const override
    {
        return 1;
    }
    virtual inline void next(size_t k, double x, double &T, double &Tp1) const override
    {
        Tp1 = boost::math::chebyshev_next(x, T, Tp1);
    }
    virtual inline void nextDerivative(size_t k, double x, double &T, double &_, double &dTp1) const override
    {
        dTp1 = (k == 0) ? 4*x : (k+2)*(2*T + dTp1/k);
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
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) const override
    {
        arma::vec c(k/2+1, arma::fill::none), f(k+1, arma::fill::zeros), df(k, arma::fill::zeros);
        for (size_t s = 0; s <= k/2; s++)
            c(s) = pow(-1, s)*pow(2, k-2*s-1)*k/(k-s)*bi(k-s, s);
        if (k%2 == 0)
            c(k/2) = pow(-1, k/2);
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
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t _, size_t k, arma::umat &bi, arma::mat &bi_05) const override
    {
        arma::vec c(k+1, arma::fill::zeros), d(k+1, arma::fill::zeros);
        c(0) = 1;
        return std::make_pair(c, d);
    }
    virtual inline double weightFunction(double _) const override
    {
        return 1;
    }
    virtual inline double weightFunctionDerivative(double _) const override
    {
        return 0;
    }
    friend class Wing;
};