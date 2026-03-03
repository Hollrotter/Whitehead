#include "B_Spline.hpp"

int main()
{
    switch (1)
    {
        case 0: // B-Spline
        {
            std::ofstream file("plot/Data/Spline/B_Spline");
            size_t M = 10;
            size_t order = 3;
            arma::vec X = arma::linspace(0, M-1, M);
            B_Spline B(X, order);

            size_t n = 200;
            arma::vec x = arma::linspace(0, M-1, n);
            arma::mat z = B(x);
            for (size_t i = 0; i < n; i++)
            {
                file << x(i) << ' ';
                for (size_t j = 0; j < M-1+order; j++)
                    file << z(j, i) << ' ';
                file << '\n';
            }
            break;
        }
        case 1:
        {
            size_t order = 3;
            size_t M = 7;
            arma::vec X = arma::linspace(0, M-1, M);
            B_Spline B1(X, order);
            B_Spline B2(X, order+1);
            for (size_t j = 0; j < M+order-1; j++) // Loop over Splines
            {
                for (size_t i = 0; i < M-1; i++) // Loop over Interval
                    std::cout << B1.integrate(i, j, order-1) << ' ';
                std::cout << " = " << B1.integrate(0.0, M-1.0, j, order-1) << '\n';
            }
            std::cout << "\n\n";
            for (size_t j = 1; j < M+order; j++) // Loop over Splines
            {
                for (size_t i = 0; i < M-1; i++) // Loop over Interval
                    std::cout << B2.integrate(i, j, order-1) << ' ';
                std::cout << " = " << B2.integrate(0.0, M-1.0, j, order-1) << '\n';
            }
            break;
        }
    }
}