#include "Lagrange.hpp"

int main()
{
    switch (2)
    {
        case 0: // Interpolation 1D
        {
            size_t n = 17; // Number of interpolation nodes
            size_t N = 100; // Number of evaluation nodes
            arma::vec x = Chebyshev::gaussLobatto(n);
            arma::vec X = arma::linspace(-1, 1, N);
            arma::mat T = Lagrange::interpolationMatrix(x, X);
            arma::vec f = 10*sin(15*x)/(15*x); // sinc
            arma::vec F = T*f;
            std::ofstream file_1("plot/Data/Lagrange/interpolation1D_1");
            for (size_t i = 0; i < n; i++)
                file_1 << x(i) << ' ' << f(i) << '\n';
            file_1.close();
            std::ofstream file_2("plot/Data/Lagrange/interpolation1D_2");
            for (size_t i = 0; i < N; i++)
                file_2 << X(i) << ' ' << F(i) << '\n';
            file_2.close();
            break;
        }
        case 1: // Interpolation 2D
        {
            size_t n = 10;
            size_t m = 30;
            arma::vec x = Chebyshev::gaussLobatto(n);
            arma::vec y = Chebyshev::gaussLobatto(m);
            arma::mat z = sin(arma::datum::pi*x)*sin(arma::datum::pi*y).t();
            arma::vec xi  = arma::linspace(-1, 1, 100);
            arma::vec eta = arma::linspace(-1, 1, 90);
            arma::mat Tx = Lagrange::interpolationMatrix(x, xi);
            arma::mat Ty = Lagrange::interpolationMatrix(y, eta);
            arma::mat Z = Lagrange::interpolation2D(Tx, Ty, z, xi, eta);

            std::ofstream file("plot/Data/Lagrange/interpolation2D");
            for (size_t i = 0; i < xi.size(); i++, file<<'\n')
                for (size_t j = 0; j < eta.size(); j++, file<<'\n')
                    file << xi(i) << ' ' << eta(j) << ' ' << Z(i, j);
            file.close();
            break;
        }
        case 2: // Interpolation curvilinear
        {
            size_t n = 20;
            size_t m = 30;
            size_t N = 5*n;
            size_t M = 5*m;
            double b = 10;
            double c = 1.39;

            arma::vec x1 = Chebyshev::gaussLobatto(n);
            arma::vec x2 = Chebyshev::gaussLobatto(m);

            arma::vec x1S = c/2*(1 + Chebyshev::gaussLobatto(n));
            arma::vec x2S = 0.9*b/4*(1 + Chebyshev::gaussLobatto(m));

            arma::mat yS = arma::ones(x1S.size(), 1) * x2S.t();
            arma::mat xS = (x1S - c/4)*sqrt(1 - pow(x2S/(b/2), 2)).t();

            arma::vec y1S = c*arma::linspace(0, 1, N);
            arma::vec y2S = 0.9*b/2*arma::linspace(0, 1, M);

            arma::mat yC = arma::ones(y1S.size(), 1)*y2S.t();
            arma::mat xC = (y1S - c/4)*sqrt(1 - pow(y2S/(b/2), 2)).t();

            arma::mat x3 = exp(-x1)*sin(arma::datum::pi*x2).t();

            arma::mat  xS_hat  = Chebyshev::DiscreteChebyshevTransform(x1, x2, xS);
            arma::mat  yS_hat  = Chebyshev::DiscreteChebyshevTransform(x1, x2, yS);
            arma::cube x_w_bar = Chebyshev::ChebyshevDerivativeCoefficients(xS_hat);
            arma::cube y_w_bar = Chebyshev::ChebyshevDerivativeCoefficients(yS_hat);

            arma::mat XC = repelem(arma::linspace(-1, 1, y1S.size()),     1, y2S.size());
            arma::mat YC = repelem(arma::linspace(-1, 1, y2S.size()).t(), y1S.size(), 1);

            #pragma omp parallel for
            for (size_t i = 0; i < y1S.size(); i++)
                for (size_t j = 0; j < y2S.size(); j++)
                    for (size_t iter = 0; iter < 50; iter++)
                    {
                        double y_w = 0;
                        double x_w = 0;
                        arma::mat J(2, 2, arma::fill::zeros);
                        arma::vec dx = {xC(i, j), yC(i, j)};
                        for (size_t p = 0; p < n; p++)
                        {
                            double Tx = boost::math::chebyshev_t(p,-XC(i, j));
                            for (size_t q = 0; q < m; q++)
                            {
                                double Txy = Tx*boost::math::chebyshev_t(q,-YC(i, j));
                                x_w += xS_hat(p, q)*Txy;
                                y_w += yS_hat(p, q)*Txy;
                                J += Txy*arma::mat{{x_w_bar(p, q, 0), x_w_bar(p, q, 1)},
                                                   {y_w_bar(p, q, 0), y_w_bar(p, q, 1)}};

                                dx -= Txy*arma::vec{xS_hat(p, q), yS_hat(p, q)};
                            }
                        }
                        dx = solve(J, dx);
                        XC(i, j) -= 0.2*dx(0);
                        YC(i, j) -= 0.2*dx(1);
                        if (norm(dx) < 1e-6)
                            break;
                    }
            arma::mat T = Lagrange::interpolationMatrix(x1, x2, XC, YC);

            arma::mat y3 = arma::reshape(T*vectorise(x3), y1S.size(), y2S.size());

            std::ofstream file("plot/Data/Lagrange/interpolationCurvilinear");
            for (size_t i = 0; i < y1S.size(); i++, file<<'\n')
                for (size_t j = 0; j < y2S.size(); j++, file<<'\n')
                    file << xC(i, j) << ' ' << yC(i, j) << ' ' << y3(i, j);
            file.close();
            break;
        }
        case 3: // Parametrization by arc length
        {
            // Cardioid
            double a = 1;
            size_t n = 100;
            arma::vec t = Chebyshev::gaussLobatto(n);
            arma::vec t_0 = arma::datum::pi/2*(t*1.01);
            arma::vec x = 2*a*(1 - cos(t_0))%cos(t_0);
            arma::vec y = 2*a*(1 - cos(t_0))%sin(t_0);

            Lagrange::CurveInterpolant chi(x, y);
            std::ofstream file("plot/Data/Lagrange/Cardioid");
            for (size_t i = 0; i < n; i++)
            {
                auto [x_chi, y_chi] = chi.evaluate(chi.getNodes().at(i));
                file << x(i) << ' ' << y(i) << ' ' << x_chi << ' ' << y_chi << '\n';
            }
            file.close();
            break;
        }
        case 4: // Transfinite Interpolation
        {
            size_t n = 10;
            size_t m = 5;

            Point p1(-1,-1);
            Point p2( 1,-1);
            Point p3( 1, 1);
            Point p4(-1, 1);

            Lagrange::CurveInterpolant chi1(p1, p2, n);
            Lagrange::CurveInterpolant chi2(p2, p3, m);
            Lagrange::CurveInterpolant chi3(p3, p4, n);
            Lagrange::CurveInterpolant chi4(p4, p1, m);

            auto [x, y] = Lagrange::TransfiniteQuadMap({&chi1, &chi2, &chi3, &chi4});
            auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics({&chi1, &chi2, &chi3, &chi4});

            std::ofstream file("plot/Data/Lagrange/TransfiniteInterpolation");
            for (size_t i = 0; i < n; i++, file<<'\n')
                for (size_t j = 0; j < m; j++, file<<'\n')
                    file << x(i, j) << ' ' << y(i, j) << ' ' <<
                    dxdx1(i, j) << ' ' << dxdx2(i, j) << ' ' << dydx1(i, j) << ' ' << dydx2(i, j);
            file.close();
            break;
        }
    }
}