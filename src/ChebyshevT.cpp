#include "ChebyshevT.hpp"

ChebyshevT::ChebyshevT(size_t _n, size_t _) : BasisFunction{_n, 0}
{
    basis = Basis::T;
    xi = Chebyshev::gauss(n);
    xg.zeros(n);
    wg.zeros(n);
    dxg.zeros(n);
    dwg.zeros(n);
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        fastgl::QuadPair gl = fastgl::GLPair(n, i+1);
        xg(i) =-gl.x();
        wg(i) = gl.weight;
    }
}

std::pair<arma::field<arma::mat>, arma::field<arma::mat>> ChebyshevT::powerSeries() const
{
    arma::field<arma::mat> f(n), df(n);
    #pragma omp parallel for
    for (size_t k = 0; k < n; k++)
    {
        f(k).zeros(n, k+1);
        df(k).zeros(n, k);
        arma::vec c(k/2+1, arma::fill::none);
        for (size_t s = 0; s <= k/2; s++)
            c(s) = pow(-1, s)*pow(2, k-2*s-1)*k/(k-s)*bi[k-s][s];
        if (k%2 == 0)
            c(k/2) = pow(-1, k/2);
        for (size_t s = 0; s <= k/2; s++)
        {
            arma::mat F(n, k+1, arma::fill::zeros);
            for (size_t t = 0; t <= k-2*s; t++)
                F.col(t) += bi[k-2*s][t] * pow(xi, k-2*s-t);
            f(k) += c(s) * F;
        }
        if (k > 0)
            for (size_t s = 0; s <= (k-1)/2; s++)
            {
                arma::mat F(n, k, arma::fill::zeros);
                for (size_t t = 0; t <= k-2*s-1; t++)
                    F.col(t) += bi[k-2*s-1][t]*pow(xi, k-2*s-1-t);
                df(k) += (k-2*s)*c(s) * F;
            }
    }
    return std::make_pair(f, df);
}

std::pair<arma::mat, arma::mat> ChebyshevT::powerSeriesWeight() const
{
    arma::mat c(n, m+1, arma::fill::zeros), d(n, m+1, arma::fill::zeros);
    c.col(0).ones();
    return std::make_pair(c, d);
}