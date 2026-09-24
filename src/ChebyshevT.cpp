#include "ChebyshevT.hpp"

ChebyshevT::ChebyshevT(size_t _n, size_t _) : BasisFunction{_n, 0}
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

std::pair<arma::vec, arma::vec> ChebyshevT::powerSeries(size_t k, double x, arma::umat &bi) const
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

std::pair<arma::vec, arma::vec> ChebyshevT::powerSeriesWeight(size_t _, size_t k, arma::umat &bi, arma::mat &bi_05) const
{
    arma::vec c(k+1, arma::fill::zeros), d(k+1, arma::fill::zeros);
    c(0) = 1;
    return std::make_pair(c, d);
}