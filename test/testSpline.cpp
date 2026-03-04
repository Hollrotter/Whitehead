#include "Splinefit.hpp"

int main()
{
    switch (3)
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
            file.close();
            break;
        }
        case 1: // Integral
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
        case 2: // Splinefit for derivative
        {
            int n = 25;
            arma::vec x = arma::linspace(0, 1, n);
            arma::vec z = cos(arma::datum::pi*x);
            arma::vec dz = -arma::datum::pi*sin(arma::datum::pi*x);
            Splinefit s(x, z, 5);
            arma::vec dzdx = s.diff(x);
            std::ofstream file("plot/Data/Spline/Splinefit");
            for (int i = 0; i < n; i++)
                file << x(i) << ' ' << dz(i) << ' ' << dzdx(i) << '\n';
            break;
        }
        case 3: // Splinefit 2D
        {
            arma::vec x1 = arma::linspace(0, 1, 20);
            arma::vec x2 = arma::linspace(0, 1, 40);

            // arma::mat z = x1*arma::ones(1, x2.size()) + arma::ones(x1.size(), 1)*x2.t();
            arma::mat z = sin(arma::datum::pi*x1)*cos(arma::datum::pi*x2).t();
            Splinefit s(x1, x2, z, 3, 7);

            // arma::mat dzdx1 = arma::ones(x1.size(), x2.size());
            // arma::mat dzdx2 = arma::ones(x1.size(), x2.size());
            arma::mat dzdx1 = arma::datum::pi*cos(arma::datum::pi*x1)*cos(arma::datum::pi*x2).t();
            arma::mat dzdx2 =-arma::datum::pi*sin(arma::datum::pi*x1)*sin(arma::datum::pi*x2).t();

            arma::vec X1 = arma::linspace(0, 1, 60);
            arma::vec X2 = arma::linspace(0, 1, 60);
            arma::mat Z = s(X1, X2);

            arma::cube gradZ = s.diff(X1, X2);
            arma::mat dZdX1 = gradZ.slice(0);
            arma::mat dZdX2 = gradZ.slice(1);

            std::ofstream file("plot/Data/Spline/Splinefit2D");
            for (int i = 0; i < x1.size(); i++)
                for (int j = 0; j < x2.size(); j++, file<<'\n')
                    file << x1(i) << ' ' << x2(j) << ' ' << z(i, j) << ' ' << dzdx1(i, j) << ' ' << dzdx2(i, j);
            file.close();

            std::ofstream file_z("plot/Data/Spline/z");
            for (int i = 0; i < X1.size(); i++, file_z<<'\n')
                for (int j = 0; j < X2.size(); j++, file_z<<'\n')
                    file_z << X1(i) << ' ' << X2(j) << ' ' << Z(i, j) << ' ' << dZdX1(i, j) << ' ' << dZdX2(i, j);
            file_z.close();
            break;
        }
    }
}