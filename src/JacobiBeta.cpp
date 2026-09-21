#include "JacobiBeta.hpp"

std::pair<arma::vec, arma::vec> JacobiBeta::powerSeries(size_t k, double x, arma::umat &bi) const
{
    arma::vec c(k+1), f(k+1), df(k);
    for (size_t pp = 0; pp <= k; pp++)
    {
        double g = bi(k, pp)*tgamma(1.5+k+pp)/tgamma(pp+1)/pow(-2., pp);
        for (size_t kk = 0; kk <= pp; kk++)
            c(kk) += g * bi(pp, kk)*pow(-1, kk);
    }
    c *= tgamma(k+1)/factorial(k)/tgamma(k+1.5);
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

std::pair<arma::vec, arma::vec> JacobiBeta::powerSeriesWeight(size_t i, size_t k, arma::umat &bi, arma::mat &bi_05) const
{
    arma::vec c(k+1), d(k+1);
    for (size_t ii = 0; ii <= k; ii++)
    {
        c(ii) = bi_05(1, ii) / pow(1+xi(i), ii-0.5);
        d(ii) = bi_05(0, ii) / pow(1+xi(i), ii+0.5)/2;
    }
    return std::make_pair(c, d);
}