#include "JacobiAlpha.hpp"

JacobiAlpha::JacobiAlpha(size_t _n, size_t _m) : BasisFunction{_n, _m}
{
    basis = Basis::PA;
    std::tie( xg,  wg) = convert<NMAX>(n, GAUJACA[n]);
    std::tie(dxg, dwg) = convert<NMAX>(n, GAUJACA2[n]);
    xi = xg;
};

std::pair<arma::field<arma::mat>, arma::field<arma::mat>> JacobiAlpha::powerSeries() const
{
    arma::field<arma::mat> f(n), df(n);
    #pragma omp parallel for
    for (size_t k = 0; k < n; k++)
    {
        f(k).zeros(n, k+1);
        df(k).zeros(n, k);
        arma::vec c(k+1);
        for (size_t pp = 0; pp <= k; pp++)
        {
            double g = bi[k][pp]*exp(LGAMMA[2+2*k+2*pp] - LGAMMA[2+2*pp] - LGAMMA[2*k+1])/pow(-2., pp);
            for (size_t kk = 0; kk <= pp; kk++)
                c(kk) += g * bi[pp][kk]*pow(-1, kk);
        }
        for (size_t s = 0; s <= k; s++)
        {
            arma::mat F(n, k+1);
            for (size_t t = 0; t <= s; t++)
                F.col(t) += bi[s][t] * pow(xi, s-t);
            f(k) += c(s) * F;
        }
        if (k > 0)
            for (size_t s = 0; s < k; s++)
            {
                arma::mat F(n, k);
                for (size_t t = 0; t <= s; t++)
                    F.col(t) += bi[s][t] * pow(xi, s-t);
                df(k) += (s+1) * c(s+1) * F;
            }
    }
    return std::make_pair(f, df);
}

std::pair<arma::mat, arma::mat> JacobiAlpha::powerSeriesWeight() const
{
    arma::mat c(n, m+1), d(n, m+1);
    #pragma omp parallel for
    for (size_t k = 0; k <= m; k++)
    {
        c.col(k) = bi_05[1][k] * pow(-1, k) / pow(1-xi, k-0.5);
        d.col(k) =-bi_05[0][k] * pow(-1, k) / pow(1-xi, k+0.5)/2;
    }
    return std::make_pair(c, d);
}