#include "misc.hpp"
#include <chrono>

int main()
{
    size_t n = 20;
    size_t m = 20;
    arma::mat A(  n,   m, arma::fill::randu);
    arma::mat B(n*m, n*m, arma::fill::randu);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 100; i++)
        arma::mat C1 = faceSplitting(B, A);
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 100; i++)
        arma::mat C2 = A^B;
    auto t3 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 100; i++)
        arma::mat C3 = repelem(vectorise(A), 1, n*m)%B;
    auto t4 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 100; i++)
        arma::mat C4 = diagmat(vectorise(A))*B;
    auto t5 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dt12 = t2-t1;
    std::chrono::duration<double> dt23 = t3-t2;
    std::chrono::duration<double> dt34 = t4-t3;
    std::chrono::duration<double> dt45 = t5-t4;

    std::cout << dt12.count() << '\n';
    std::cout << dt23.count() << '\n';
    std::cout << dt34.count() << '\n';
    std::cout << dt45.count() << '\n';
}