#include "ChebyshevU.hpp"

ChebyshevU::ChebyshevU(size_t _n, size_t _m) : BasisFunction{_n, _m}
{
    xi =-cos(arma::datum::pi*(arma::regspace(0, n-1)+1)/(n+1));
    xg = xi;
    wg = arma::datum::pi/(n+1)*(1-xi%xi);
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

std::pair<arma::vec, arma::vec> ChebyshevU::powerSeriesWeight(size_t i, arma::mat &bi_05) const
{
    arma::vec c(m+1), d(m+1);
    for (size_t ii = 0; ii <= m; ii++)
        for (size_t jj = 0; jj <= ii/2; jj++)
            c(ii) += pow(-1, ii-jj)*bi_05(1, ii-jj)*bi[ii-jj][jj]*pow(2*xi(i), ii-2*jj)/pow(1-pow(xi(i), 2), ii-jj);
    c *= sqrt(1-pow(xi(i), 2));
    for (size_t ii = 0; ii < m; ii++)
        for (size_t jj = 0; jj <= ii/2; jj++)
        {
            double con = pow(-1, ii-jj)*bi[ii-jj][jj]*pow(2*xi(i), ii-2*jj)/pow(1-pow(xi(i), 2), ii-jj);
            d(ii)   += con * bi_05(0, ii-jj) * xi(i);
            d(ii+1) += con * bi_05(0, ii-jj);
        }
    d /=-sqrt(1-pow(xi(i), 2));
    return std::make_pair(c, d);
}