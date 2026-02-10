#include "Chebyshev.hpp"

arma::cx_vec Chebyshev::FFFTEO(const arma::vec f)
{
    const size_t N = f.size();
    const size_t N_2 = N/2;
    // w_p multiplied by imaginary unit, because it is only needed as this product!
    const arma::cx_double i(0, 1);
    const arma::cx_vec w_p = i*exp(-2*arma::datum::pi/N*i*arma::regspace(0, N_2));
    arma::cx_vec Z(N_2);
    for (size_t j = 0; j < N_2; j++)
        Z(j) = f(2*j) + i*f(2*j+1);
    Z = arma::fft(Z);
    arma::cx_vec F(N_2+1);
    F(0)   = (Z(0).real() + Z(0).imag())/N;
    F(N_2) = (Z(0).real() - Z(0).imag())/N;
    for (size_t k = 1; k < N_2; k++)
        F(k) = (Z(k) + conj(Z(N_2-k)) - w_p(k)*(Z(k) - conj(Z(N_2-k))))/(2.0*N);
    return F;
}

std::pair<arma::vec, arma::vec> Chebyshev::forwardRealFFT(const arma::vec x)
{
    const arma::cx_vec X = FFFTEO(x);
    const arma::vec a = 2*real(X);
    arma::vec       b =-2*imag(X);
    b(0)     = 0;
    b.back() = 0;
    return std::tie(a, b);
}

arma::vec Chebyshev::fastCosineTransform(const arma::vec f, const std::string s)
{
    const size_t N = f.size()-1;
    const arma::vec C = cos(arma::regspace(0, N)*arma::datum::pi/N);
    const arma::vec S = sin(arma::regspace(0, N)*arma::datum::pi/N);
    arma::vec e(N);
    for (size_t j = 0; j < N; j++)
        e(j) = (f(j) + f(N-j))/2 - S(j)*(f(j) - f(N-j));
    auto [a_bar, b_bar] = forwardRealFFT(e);
    const size_t N_2 = N/2;
    arma::vec a(N+1);
    for (size_t k = 0; k < N_2+1; k++)
        a(2*k) = a_bar(k);
    a(1) = f(0) - f.back();
    for (size_t j = 1; j < N; j++)
        a(1) += 2*C(j)*f(j);
    a(1) /= N;
    for (size_t k = 1; k < N_2; k++)
        a(2*k+1) = b_bar(k) + a(2*k-1);
    if (s == "BACKWARD")
        a *= N_2;
    return a;
}

arma::vec Chebyshev::fastChebyshevTransform(const arma::vec f, const std::string s)
{
    arma::vec g = f;
    
    if (s == "BACKWARD")
    {
        g(0)     *= 2;
        g.back() *= 2;
    }

    arma::vec a = fastCosineTransform(g, s);

    if (s == "FORWARD")
    {
        a(0)     /= 2;
        a.back() /= 2;
    }

    return a;
}

arma::mat Chebyshev::fastChebyshevTransform(const arma::mat f, const std::string s)
{
    arma::mat f_hat = f;

    #pragma omp parallel for
    for (size_t j = 0; j < f.n_cols; j++)
        f_hat.col(j) = fastChebyshevTransform(arma::vec(f_hat.col(j)), s);

    f_hat = f_hat.t();

    #pragma omp parallel for
    for (size_t i = 0; i < f.n_rows; i++)
        f_hat.col(i) = fastChebyshevTransform(arma::vec(f_hat.col(i)), s);

    return f_hat.t();
}