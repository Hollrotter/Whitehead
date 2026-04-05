#include "interpolation.hpp"

arma::mat ChebyshevInterpolation(arma::vec x1, arma::vec x2, arma::mat xs, arma::mat ys, arma::mat xss, arma::mat yss)
{
    size_t nx = x1.size();
    size_t ny = x2.size();
    arma::mat x_m_hat, y_m_hat;
    if (nx%2 == 0 && ny%2 == 0)
    {
        x_m_hat = Chebyshev::DiscreteChebyshevTransform(x1, x2, xs);
        y_m_hat = Chebyshev::DiscreteChebyshevTransform(x1, x2, ys);
    }
    else
    {
        x_m_hat = Chebyshev::fastChebyshevTransform(xs, "FORWARD");
        y_m_hat = Chebyshev::fastChebyshevTransform(ys, "FORWARD");
    }
    arma::cube x_m_bar = Chebyshev::ChebyshevDerivativeCoefficients(x_m_hat);
    arma::cube y_m_bar = Chebyshev::ChebyshevDerivativeCoefficients(y_m_hat);

    arma::mat xC = repelem(Chebyshev::gaussLobatto(nx),     1, ny);
    arma::mat yC = repelem(Chebyshev::gaussLobatto(ny).t(), nx, 1);
    
    #pragma omp parallel for
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
            for (size_t iter = 0; iter < 1000; iter++)
            {
                arma::mat J(2, 2, arma::fill::zeros);
                arma::vec b = {xss(i, j), yss(i, j)};
                for (size_t p = 0; p < nx; p++)
                {
                    double Tx = boost::math::chebyshev_t(p,-xC(i, j));
                    for (size_t q = 0; q < ny; q++)
                    {
                        double Txy = Tx*boost::math::chebyshev_t(q,-yC(i, j));
                        J += Txy*arma::mat{{x_m_bar(p, q, 0), x_m_bar(p, q, 1)},
                                           {y_m_bar(p, q, 0), y_m_bar(p, q, 1)}};
                        b -= Txy*arma::vec{x_m_hat(p, q), y_m_hat(p, q)};
                    }
                }
                b = solve(J, b);
                xC(i, j) -= b(0);
                yC(i, j) -= b(1);
                if (norm(b) < 1e-10)
                    break;
            }
    return Lagrange::interpolationMatrix(x1, x2, xC, yC);
}

arma::mat LagrangeInterpolation(arma::vec x1, arma::vec x2, arma::mat xs, arma::mat ys, arma::mat xss, arma::mat yss)
{
    size_t nx = x1.size();
    size_t ny = x2.size();

    arma::mat Tx1 = Lagrange::interpolationMatrix(x1, xs);
    arma::mat Ty1 = Lagrange::interpolationMatrix(x1, ys);
    arma::mat Tx2 = Lagrange::interpolationMatrix(x2, xs);
    arma::mat Ty2 = Lagrange::interpolationMatrix(x2, ys);

    arma::mat Tx = kron(Tx1, Tx2);
    arma::mat Ty = kron(Ty1, Ty2);

    arma::mat D1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
    arma::mat D2 = Chebyshev::derivativeMatrix(x2, Derivative::first);

    arma::mat xC = repelem(Chebyshev::gaussLobatto(nx),     1, ny);
    arma::mat yC = repelem(Chebyshev::gaussLobatto(ny).t(), nx, 1);
    
    #pragma omp parallel for
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
            for (size_t iter = 0; iter < 1000; iter++)
            {
                double dxCd1 = dot(D1.row(i), xC.col(j));
                double dyCd1 = dot(D1.row(i), yC.col(j));
                double dxCd2 = dot(xC.row(i), D2.row(j));
                double dyCd2 = dot(yC.row(i), D2.row(j));
                arma::mat J(2, 2, arma::fill::zeros);
                arma::vec b = {xss(i, j), yss(i, j)};
                for (size_t p = 0; p < nx; p++)
                {
                    for (size_t q = 0; q < ny; q++)
                    {
                        J += arma::mat{{dxCd1, dxCd2},
                                       {dyCd1, dyCd2}};
                        b -= arma::vec{dot(Tx.row(q+ny*p), vectorise(xC)), dot(Ty.row(1+ny*p), vectorise(yC))};
                    }
                }
                b = solve(J, b);
                xC(i, j) -= b(0);
                yC(i, j) -= b(1);
                if (norm(b) < 1e-10)
                    break;
            }
    return Lagrange::interpolationMatrix(x1, x2, xC, yC);
}