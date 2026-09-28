#include "JacobiAlpha.hpp"

JacobiAlpha::JacobiAlpha(size_t _n, size_t _m) : BasisFunction{_n, _m}
{
    basis = Basis::PA;
    std::tie(xg, wg) = gaujac(n, 0.5, 0.0);
    xi = xg;
    std::tie(dxg, dwg) = gaujac(n, -0.5, 0.0);
};

std::pair<arma::vec, arma::vec> JacobiAlpha::powerSeries(size_t k, double x) const
{
    arma::vec c(k+1), f(k+1), df(k);
    for (size_t pp = 0; pp <= k; pp++)
    {
        double g = bi[k][pp]*tgamma(1.5+k+pp)/tgamma(1.5+pp)/pow(-2., pp);
        for (size_t kk = 0; kk <= pp; kk++)
            c(kk) += g * bi[pp][kk]*pow(-1, kk);
    }
    c /= factorial(k);
    for (size_t s = 0; s <= k; s++)
    {
        arma::vec F(k+1);
        for (size_t t = 0; t <= s; t++)
            F(t) += bi[s][t] * pow(x, s-t);
        f += c(s) * F;
    }
    if (k > 0)
        for (size_t s = 0; s < k; s++)
        {
            arma::vec F(k);
            for (size_t t = 0; t <= s; t++)
                F(t) += bi[s][t] * pow(x, s-t);
            df += (s+1) * c(s+1) * F;
        }
    return std::make_pair(f, df);
}

std::pair<arma::vec, arma::vec> JacobiAlpha::powerSeriesWeight(size_t i, arma::mat &bi_05) const
{
    arma::vec c(m+1), d(m+1);
    for (size_t k = 0; k <= m; k++)
    {
        c(k) = bi_05(1, k) * pow(-1, k) / pow(1-xi(i), k-0.5);
        d(k) =-bi_05(0, k) * pow(-1, k) / pow(1-xi(i), k+0.5)/2;
    }
    return std::make_pair(c, d);
}