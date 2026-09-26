#pragma once
#include "BasisFunction.hpp"

class JacobiAlpha : public BasisFunction
{
public:
    JacobiAlpha(size_t _n, size_t _m);
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
        return (1 + 5*x)/4;
    }
    virtual inline arma::vec linear(arma::vec x) const override
    {
        return (1 + 5*x)/4;
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
        return 1.25;
    }
    virtual inline arma::vec linearDerivative(size_t k) const override
    {
        return 1.25*arma::ones(k);
    }
    virtual inline void next(size_t k, double x, double &P, double &Pp1) const override
    {
        double gamma = 2*k + 3.5;
        Pp1 = (gamma*(0.25 + (pow(gamma, 2) - 1)*x)*P - 2*(k+1.5)*(k+1)*(gamma + 1)*Pp1)
            / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline void next(size_t k, arma::vec x, arma::mat &P) const override
    {
        double gamma = 2*k + 3.5;
        P.col(k+1) = (gamma*(0.25 + (pow(gamma, 2) - 1)*x)%P.col(k) - 2*(k+1.5)*(k+1)*(gamma + 1)*P.col(k-1))
                   / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline void nextDerivative(size_t k, double x, double &P, double &dP, double &dPp1) const override
    {
        double gamma = 2*k + 3.5;
        dPp1 = (gamma*(pow(gamma, 2) - 1)*P + gamma*(0.25 + (pow(gamma, 2)-1)*x)*dP - 2*(k+1.5)*(k+1)*(gamma+1)*dPp1)
             / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline void nextDerivative(size_t k, arma::vec x, arma::mat &P, arma::mat &dP) const override
    {
        double gamma = 2*k + 3.5;
        dP.col(k+1) = (gamma*(pow(gamma, 2) - 1)*P.col(k) + gamma*(0.25 + (pow(gamma, 2)-1)*x)%dP.col(k) - 2*(k+1.5)*(k+1)*(gamma+1)*dP.col(k-1))
                    / (2*(k+2)*(gamma-k-1)*(gamma-1));
    }
    virtual inline double left(size_t k) const override
    {
        return k%2==0 ? 1 : -1;
    }
    virtual inline double right(size_t k) const override
    {
        return tgamma(k+1.5)/factorial(k)/tgamma(1.5);
    }
    virtual inline double leftDerivative(size_t k) const override
    {
        return k<1 ? 0 : k%2==0 ? -(k+1.5)*tgamma(k+1)/2/factorial(k-1)/tgamma(2) : (k+1.5)*tgamma(k+1)/2/factorial(k-1)/tgamma(2);
    }
    virtual inline double rightDerivative(size_t k) const override
    {
        return k<1 ? 0 : (k+1.5)*tgamma(k+1.5)/2/factorial(k-1)/tgamma(2.5);
    }
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x) const override;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t i, arma::mat &bi_05) const override;
    virtual inline double weightFunction(double x) const override
    {
        return sqrt(1 - x);
    }
    virtual inline arma::vec weightFunction(arma::vec x) const override
    {
        return sqrt(1 - x);
    }
    virtual inline double weightFunctionDerivative(double x) const override
    {
        return -0.5/sqrt(1 - x);
    }
    virtual inline arma::vec weightFunctionDerivative(arma::vec x) const override
    {
        return -0.5/sqrt(1 - x);
    }
    friend class Wing;
private:
    size_t factorial(size_t k) const
    {
        size_t f = 1;
        for (size_t i = 2; i <= k; i++)
            f *= i;
        return f;
    }
};