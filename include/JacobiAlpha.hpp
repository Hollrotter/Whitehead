#pragma once
#include "BasisFunction.hpp"

class JacobiAlpha : public BasisFunction
{
public:
    JacobiAlpha(size_t _n) : BasisFunction{_n} {};
    virtual inline constexpr double constant() override
    {
        return 1;
    }
    virtual inline double linear(double x) override
    {
        return (1 + 5*x)/4;
    }
    virtual inline constexpr double constantDerivative() override
    {
        return 0;
    }
    virtual inline constexpr double linearDerivative() override
    {
        return 1.25;
    }
    virtual inline void next(size_t k, double x, double &P, double &Pp1) override
    {
        double gamma = 2*k + 3.5;
        Pp1 = (gamma*(0.25 + (pow(gamma, 2) - 1)*x)*P - 2*(k+1.5)*(k+1)*(gamma + 1)*Pp1)
            / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline void nextDerivative(size_t k, double x, double &P, double &dP, double &dPp1) override
    {
        double gamma = 2*k + 3.5;
        dPp1 = (gamma*(pow(gamma, 2) - 1)*P + gamma*(0.25 + (pow(gamma, 2)-1)*x)*dP - 2*(k+1.5)*(k+1)*(gamma+1)*dPp1)
             / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline double left(size_t k) override
    {
        return k%2==0 ? 1 : -1;
    }
    virtual inline double right(size_t k) override
    {
        return tgamma(k+1.5)/factorial(k)/tgamma(1.5);
    }
    virtual inline double leftDerivative(size_t k) override
    {
        return k<1 ? 0 : k%2==0 ? -(k+1.5)*tgamma(k+1)/2/factorial(k-1)/tgamma(2) : (k+1.5)*tgamma(k+1)/2/factorial(k-1)/tgamma(2);
    }
    virtual inline double rightDerivative(size_t k) override
    {
        return k<1 ? 0 : (k+1.5)*tgamma(k+1.5)/2/factorial(k-1)/tgamma(2.5);
    }
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) override
    {
        arma::vec c(k+1), f(k+1), df(k);
        for (size_t pp = 0; pp <= k; pp++)
        {
            double g = bi(k, pp)*tgamma(1.5+k+pp)/tgamma(1.5+pp)/pow(-2., pp);
            for (size_t kk = 0; kk <= pp; kk++)
                c(kk) += g * bi(pp, kk)*pow(-1, kk);
        }
        c *= tgamma(1.5+k)/factorial(k)/tgamma(1.5+k);
        for (size_t s = 0; s <= k; s++)
        {
            arma::vec F(k+1);
            for (size_t t = 0; t <= s; t++)
                F(t) += bi(s, t) * pow(x, s-t);
            f += c(s) * F;
        }
        if (k > 0)
            for (size_t s = 0; s < k; s++)
            {
                arma::vec F(k);
                for (size_t t = 0; t <= s; t++)
                    F(t) += bi(s, t) * pow(x, s-t);
                df += (s+1) * c(s+1) * F;
            }
        return std::make_pair(f, df);
    }
    friend class Wing;
private:
    size_t factorial(size_t k)
    {
        size_t f = 1;
        for (size_t i = 2; i <= k; i++)
            f *= i;
        return f;
    }
};