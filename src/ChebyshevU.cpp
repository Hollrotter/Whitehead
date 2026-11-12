#include "ChebyshevU.hpp"

ChebyshevU::ChebyshevU(size_t _n, size_t _m) : BasisFunction{_n, _m}
{
    basis = Basis::U;
    xi =-cos(arma::datum::pi*(arma::regspace(0, n-1)+1)/(n+1));
    xg = xi;
    wg = arma::datum::pi/(n+1)*(1-xi%xi);
    std::tie(dxg, dwg) = convert<NMAX>(n, GAUJACU2[n]);
}

std::pair<arma::vec, arma::vec> ChebyshevU::powerSeries(size_t k, double x) const
{
    arma::vec c(k/2+1, arma::fill::none), f(k+1, arma::fill::zeros), df(k, arma::fill::zeros);
    for (size_t s = 0; s <= k/2; s++)
        c(s) = pow(-1, s)*pow(2, k-2*s)*bi[k-s][s];
    for (size_t s = 0; s <= k/2; s++)
    {
        arma::vec F(k+1, arma::fill::zeros);
        for (size_t t = 0; t <= k-2*s; t++)
            F(t) += bi[k-2*s][t] * pow(x, k-2*s-t);
        f += c(s) * F;
    }
    if (k > 0)
        for (size_t s = 0; s <= (k-1)/2; s++)
        {
            arma::vec F(k, arma::fill::zeros);
            for (size_t t = 0; t <= k-2*s-1; t++)
                F(t) += bi[k-2*s-1][t]*pow(x, k-2*s-1-t);
            df += (k-2*s)*c(s) * F;
        }
    return std::make_pair(f, df);
}

std::pair<arma::mat, arma::mat> ChebyshevU::powerSeriesWeight() const
{
    arma::mat c(n, m+1, arma::fill::zeros), d(n, m+1, arma::fill::zeros);
    for (size_t ii = 0; ii <= m; ii++)
    {
        for (size_t jj = 0; jj <= ii/2; jj++)
            c.col(ii) += pow(-1, ii-jj)*bi_05[1][ii-jj]*bi[ii-jj][jj]*pow(2*xi, ii-2*jj)/pow(1-pow(xi, 2), ii-jj);
        c.col(ii) %= sqrt(1-pow(xi, 2));
    }
    for (size_t ii = 0; ii < m; ii++)
    {
        for (size_t jj = 0; jj <= ii/2; jj++)
        {
            arma::vec con = pow(-1, ii-jj)*bi[ii-jj][jj]*pow(2*xi, ii-2*jj)/pow(1-pow(xi, 2), ii-jj);
            d.col(ii)   += bi_05[0][ii-jj] * con % xi;
            d.col(ii+1) += bi_05[0][ii-jj] * con;
        }
        d.col(ii) /=-sqrt(1-pow(xi, 2));
    }
    return std::make_pair(c, d);
}