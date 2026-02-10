#include "Chebyshev.hpp"

int main()
{
    switch (2)
    {
        case 0: // Derivative matrices
        {
            arma::vec x   = Chebyshev::gaussLobatto(100);
            arma::vec z   = exp(-pow(3*x, 2));
            arma::vec dz  = Chebyshev::derivativeMatrix(x, Derivative::first)  * z;
            arma::vec d2z = Chebyshev::derivativeMatrix(x, Derivative::second) * z;
            z.t().print();
            dz.t().print();
            d2z.t().print();
            std::ofstream file("plot/Data/Chebyshev/Chebyshev");
            for (size_t i = 0; i < x.size(); i++)
                file << x(i) << ' ' << z(i) << ' ' << dz(i) << ' ' << d2z(i) << '\n';
            file.close();
            break;
        }
        case 1: // Discrete Chebyshev Transform
        {
			size_t n = 7; // must be odd!
			arma::vec x = Chebyshev::gaussLobatto(n);
			arma::vec a = {1, 2, 3, 4, 5, 6, 7};
			arma::vec u(n, arma::fill::zeros);
			for (size_t i = 0; i < a.size(); i++)
				u += a(i)*Chebyshev::Polynomial(i,-x);
			arma::vec u_hat = Chebyshev::DiscreteChebyshevTransform(x, u);
			arma::vec u_hatFFT = Chebyshev::fastChebyshevTransform(u, "FORWARD");
			u_hat.t().print();
			u_hatFFT.t().print();
			std::cout << (u_hat - u_hatFFT).t() << '\n';
            break;
        }
        case 2: // Discrete Chebyshev Transform 2D
        {
			size_t n = 21;
			size_t m = 21;

			arma::vec x = Chebyshev::gaussLobatto(n);
			arma::vec y = Chebyshev::gaussLobatto(m);
			arma::mat u = exp(-2*x)*exp(-3*y).t();
			auto t1 = std::chrono::high_resolution_clock::now();
			arma::mat u_hat = Chebyshev::DiscreteChebyshevTransform(x, y, u);
			auto t2 = std::chrono::high_resolution_clock::now();
			arma::mat u_hatFFT = Chebyshev::fastChebyshevTransform(u, "FORWARD");
			auto t3 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> dt12 = t2 - t1;
			std::chrono::duration<double> dt23 = t3 - t2;

			std::cout << "DFT: " << dt12.count() << '\n';
			std::cout << "FFT: " << dt23.count() << '\n';
            break;
        }
        case 3: // Roots 1D
        {
			size_t n = 10;
			arma::vec x = Chebyshev::gaussLobatto(n);
			arma::vec u = sin(arma::datum::pi*x)/(arma::datum::pi*x);
			arma::vec u_hat = Chebyshev::DiscreteChebyshevTransform(x, u);
            printf("Function values:\n");
			u.t().print();
            printf("Chebyshev coefficients:\n");
			u_hat.t().print();

			arma::mat A(n-1, n-1, arma::fill::zeros);
			A(0, 1) = 1;
			for (size_t j = 1; j < n-2; j++)
			{
				A(j, j-1) = 0.5;
				A(j, j+1) = 0.5;
			}
			A.row(n-2) =-u_hat.head(n-1).t()/u_hat.back()/2;
			A(n-2, n-3) += 0.5;

			arma::cx_vec eigval;
			arma::cx_mat eigvec;

			arma::eig_gen(eigval, eigvec, A);
            printf("Zeros calculated by an eigenvalue problem:\n");
			eigval.print();
            break;
        }
    }
}