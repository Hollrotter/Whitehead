#pragma once
#include "BasisFunction.hpp"
#include "Chebyshev.hpp"
#include "fastgl.h"

class ChebyshevT : public BasisFunction
{
public:
    ChebyshevT(size_t _n, size_t _) : BasisFunction{_n, 0}
    {
        xi = Chebyshev::gauss(n);
        xg.zeros(n);
        wg.zeros(n);
        for (size_t i = 0; i < n; i++)
        {
            fastgl::QuadPair gl = fastgl::GLPair(n, i+1);
            xg(i) =-gl.x();
            wg(i) = gl.weight;
        }
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
    virtual std::pair<arma::vec, arma::vec> powerSeries(size_t k, double x, arma::umat &bi) const override;
    virtual std::pair<arma::vec, arma::vec> powerSeriesWeight(size_t _, size_t k, arma::umat &bi, arma::mat &bi_05) const override;
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