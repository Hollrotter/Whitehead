#include "Lagrange.hpp"

/**
 * @brief 
 * 
 * @param x Vector of nodes forming the Gauss-Lobatto Nodes.
 * @return arma::vec 
 */
arma::vec Lagrange::barycentricWeights(const arma::vec x)
{
    arma::vec w(x.size(), arma::fill::ones);
    for (size_t j = 1; j < x.size(); j++)
        for (size_t k = 0; k < j; k++)
        {
            w(k) *= x(k) - x(j);
            w(j) *= x(j) - x(k);
        }
    return arma::min(abs(w))/w;
}

arma::vec Lagrange::barycentricWeights(const size_t n)
{
    auto delta = [](const size_t i, const size_t n) { return (i == 0 || i == n-1) ? 0.5 : 1; };
    arma::vec w(n);
    for (size_t i = 0; i < n; i++)
        w(i) =-pow(-1, n-i)*delta(i, n);
    return w;
}

double Lagrange::interpolation(const double x, const arma::vec X, const arma::vec f, const arma::vec w)
{
    double numerator   = 0;
    double denominator = 0;
    for (size_t j = 0; j < f.size(); j++)
    {
        if (almostEqual(x, X(j)))
            return f(j);
        double t = w(j)/(x - X(j));
        numerator   += t*f(j);
        denominator += t;
    }
    return numerator/denominator;
}

arma::vec Lagrange::interpolation(const arma::vec x, const arma::vec X, const arma::vec f, const arma::vec w)
{
    arma::vec xi(x.size());
    for (size_t i = 0; i < x.size(); i++)
        xi(i) = interpolation(x(i), X, f, w);
    return xi;
}

double Lagrange::interpolantDerivative(const double x, const arma::vec X, const arma::vec f, const arma::vec w)
{
    size_t n = f.size();
    bool atNode = false;
    double denominator = 0;
    double numerator   = 0;
    double p = 0;
    size_t i =-1;
    for (size_t j = 0; j < n; j++)
    {
        if (almostEqual(x, X(j)))
        {
            atNode = true;
            p = f(j);
            denominator =-w(j);
            i = j;
        }
    }
    if (atNode)
    {
        for (size_t j = 0; j < n; j++)
            if (i != j)
                numerator += w(j)*(p - f(j))/(x - X(j));
    }
    else
    {
        denominator = 0;
        p = interpolation(x, X, f, w);
        for (size_t j = 0; j < n; j++)
        {
            double t = w(j)/(x - X(j));
            numerator   += t*(p - f(j))/(x - X(j));
            denominator += t;
        }
    }
    return numerator/denominator;
}

arma::vec Lagrange::interpolantDerivative(const arma::vec x, const arma::vec X, const arma::vec f, const arma::vec w)
{
    arma::vec dx(x.size());
    for (size_t i = 0; i < x.size(); i++)
        dx(i) = interpolantDerivative(x(i), X, f, w);
    return dx;
}

/**
 * @brief 
 * 
 * @param x Vector containing the coordinates of the initial mesh.
 * @param xi Vector containing the coordinates of the target mesh.
 * @return arma::mat 
 */
arma::mat Lagrange::interpolationMatrix(const arma::vec x, const arma::vec xi)
{
    size_t n = x.size();
    arma::vec w = barycentricWeights(x);
    arma::mat T(xi.size(), n, arma::fill::zeros);
    #pragma omp parallel for
    for (size_t k = 0; k < xi.size(); k++)
    {
        bool rowHasMatch = false;
        for (size_t j = 0; j < n; j++)
            if (almostEqual(xi(k), x(j)))
            {
                rowHasMatch = true;
                T(k, j) = 1;
            }
        if (rowHasMatch == false)
        {
            double s = 0;
            for (size_t j = 0; j < n; j++)
            {
                double t = w(j)/(xi(k) - x(j));
                T(k, j) = t;
                s += t;
            }
            T.row(k) /= s;
        }
    }
    return T;
}

/**
 * @brief 
 * 
 * @param x1 Vector containing nodes of 1-coordinate of computational source grid.
 * @param x2 Vector containing nodes of 2-coordinate of computational source grid.
 * @param x Matrix containing the nodes of the 1-coordinate of the computational target grid.
 * @param y Matrix containing the nodes of the 2-coordinate of the computational target grid.
 * @return arma::mat 
 */
arma::mat Lagrange::interpolationMatrix(const arma::vec x1, const arma::vec x2, const arma::mat x, const arma::mat y)
{
    // Check if the matrices x and y have the same shape:
    if (x.n_rows != y.n_rows || x.n_cols != y.n_cols)
    {
        std::println("Interpolationmatrix: size mismatch between x and y!");
        exit(EXIT_FAILURE);
    }

    // Some array lengths
    size_t n1 = x1.size();
    size_t n2 = x2.size();
    size_t nx = x.n_rows;
    size_t ny = x.n_cols;

    // Initialization of the interpolation matrix
    arma::mat T(nx*ny, n1*n2);

    // Interpolation in x1-direction
    for (size_t n = 0; n < ny; n++)
    {
        arma::mat T1 = interpolationMatrix(x1, x.col(n));
        for (size_t i = 0; i < n2; i++)
            T.rows(n*nx, (n+1)*nx-1).cols(i*n1, (i+1)*n1-1) = T1;
    }

    // Interpolation in x2-direction
    for (size_t m = 0; m < nx; m++)
    {
        arma::mat T2 = interpolationMatrix(x2, y.row(m).t());
        for (size_t i = 0; i < ny; i++)
            for (size_t j = 0; j < n2; j++)
                T.row(i*nx+m).cols(j*n1, (j+1)*n1-1) *= T2(i, j);
    }
    
    return T;
}

/**
 * @brief 
 * 
 * @param Tx Matrix for the interpolation in x-direction.
 * @param Ty Matrix for the interpolation in y-direction.
 * @param z Matrix of function values at the given grid points.
 * @param xi Vector containing the nodes of the target mesh in x-direction.
 * @param eta Vector containing the nodes of the target mesh in y-direction.
 * @return arma::mat 
 */
arma::mat Lagrange::interpolation2D(const arma::mat Tx, const arma::mat Ty, const arma::mat z, const arma::vec xi, const arma::vec eta)
{
    arma::mat Z_bar(xi.size(), Ty.n_cols);
    for (size_t j = 0; j < Ty.n_cols; j++)
        Z_bar.col(j) = Tx*z.col(j);
    arma::mat Z(xi.size(), eta.size());
    for (size_t n = 0; n < xi.size(); n++)
        Z.row(n) = Z_bar.row(n)*Ty.t();
    return Z;
}

arma::vec Lagrange::CurveInterpolant::parametrize()
{
    size_t n = x.size();
    arma::vec s = Chebyshev::gaussLobatto(n);
    arma::vec t(n, arma::fill::zeros);
    arma::mat D = Chebyshev::derivativeMatrix(s, Derivative::first);
    arma::vec dxds = D*x;
    arma::vec dyds = D*y;

    arma::vec I = sqrt(pow(dxds, 2) + pow(dyds, 2));
    D.row(0).eye();
    I(0) = 0;
    arma::vec int_I_dx = solve(D, I);
    double L = int_I_dx(n-1);

    arma::vec w = barycentricWeights(n);
    t(0) =-1;
    double I2 = sqrt(pow(interpolantDerivative(-1, s, x, w), 2)
                   + pow(interpolantDerivative(-1, s, y, w), 2));
    double I1;
    for (size_t i = 1; i < n-1; i++)
    {
        t(i) = t(i-1) + s(i) - s(i-1);
        for (size_t q = 0; q < 1000; q++)
        {
            double dxdt = interpolantDerivative((t(i)+t(i-1))/2, s, x, w);
            double dydt = interpolantDerivative((t(i)+t(i-1))/2, s, y, w);

            double Im = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            dxdt = interpolantDerivative(t(i), s, x, w);
            dydt = interpolantDerivative(t(i), s, y, w);

            I1 = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            double f = L/2*(s(i) - s(i-1)) - (t(i) - t(i-1))*(I1 + 4*Im + I2)/6;
            t(i) += 0.1*f/I1;
            if (fabs(f) < 1e-10)
                break;
        }
        I2 = I1;
    }
    t.back() = 1;
    return t;
}

std::pair<arma::mat, arma::mat> Lagrange::TransfiniteQuadMap(const std::array<CurveInterpolant*, 4> chi)
{
    if (chi[0]->getNodes().size() != chi[2]->getNodes().size() || chi[1]->getNodes().size() != chi[3]->getNodes().size())
    {
        std::println("TransfiniteQuadMap: size mismatch!");
        exit(EXIT_FAILURE);
    }
    double sign1, sign2, sign3, sign4;
    if (almostEqual(chi[0]->evaluate(-1.0), chi[3]->evaluate(-1.0)))
    {
        sign1 = 1.;
        sign4 = 1.;
    }
    else if (almostEqual(chi[0]->evaluate(1.0), chi[3]->evaluate(-1.0)))
    {
        sign1 =-1.;
        sign4 = 1.;
    }
    else if (almostEqual(chi[0]->evaluate(1.0), chi[3]->evaluate(1.0)))
    {
        sign1 =-1.;
        sign4 =-1.;
    }
    else if (almostEqual(chi[0]->evaluate(-1.0), chi[3]->evaluate(1.0)))
    {
        sign1 = 1.;
        sign4 =-1.;
    }
    else
    {
        std::println("TransfiniteQuadMap: chi1 and chi4 mismatch!");
        exit(EXIT_FAILURE);
    }
    if (almostEqual(chi[1]->evaluate(1.0), chi[2]->evaluate(1.0)))
    {
        sign2 = 1;
        sign3 = 1;
    }
    else if (almostEqual(chi[1]->evaluate(-1.0), chi[2]->evaluate(1.0)))
    {
        sign2 =-1;
        sign3 = 1;
    }
    else if (almostEqual(chi[1]->evaluate(-1.0), chi[2]->evaluate(-1.0)))
    {
        sign2 =-1;
        sign3 =-1;
    }
    else if (almostEqual(chi[1]->evaluate(1.0), chi[2]->evaluate(-1.0)))
    {
        sign2 = 1;
        sign3 =-1;
    }
    else
    {
        std::println("TransfiniteQuadMap: chi2 and chi3 mismatch!");
        exit(EXIT_FAILURE);
    }
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    arma::mat x(x1.size(), x2.size());
    arma::mat y(x1.size(), x2.size());
    for (size_t i = 0; i < x1.size(); i++)
    {
        double x1i = x1(i);
        auto [X_1, Y_1] = chi[0]->evaluate(sign1*x1i);
        auto [X_3, Y_3] = chi[2]->evaluate(sign3*x1i);
        for (size_t j = 0; j < x2.size(); j++)
        {
            double x2j = x2(j);
            auto [X_2, Y_2] = chi[1]->evaluate(sign2*x2j);
            auto [X_4, Y_4] = chi[3]->evaluate(sign4*x2j);
            x(i, j) = ((1 - x1i)*X_4 + (1 + x1i)*X_2 + (1 - x2j)*X_1 + (1 + x2j)*X_3)/2
                    - ((1 - x1i)*((1 - x2j)*x_1 + (1 + x2j)*x_4) + (1 + x1i)*((1 - x2j)*x_2 + (1 + x2j)*x_3))/4;
            y(i, j) = ((1 - x1i)*Y_4 + (1 + x1i)*Y_2 + (1 - x2j)*Y_1 + (1 + x2j)*Y_3)/2
                    - ((1 - x1i)*((1 - x2j)*y_1 + (1 + x2j)*y_4) + (1 + x1i)*((1 - x2j)*y_2 + (1 + x2j)*y_3))/4;
        }
    }
    return std::tie(x, y);
}

std::tuple<arma::mat, arma::mat, arma::mat, arma::mat> Lagrange::TransfiniteQuadMetrics(const std::array<CurveInterpolant*, 4> chi)
{
    if (chi[0]->getNodes().size() != chi[2]->getNodes().size() || chi[1]->getNodes().size() != chi[3]->getNodes().size())
    {
        std::println("TransfiniteQuadMetric: size mismatch!");
        exit(EXIT_FAILURE);
    }
    double sign1, sign2, sign3, sign4;
    if (almostEqual(chi[0]->evaluate(-1.0), chi[3]->evaluate(-1.0)))
    {
        sign1 = 1.;
        sign4 = 1.;
    }
    else if (almostEqual(chi[0]->evaluate(1.0), chi[3]->evaluate(-1.0)))
    {
        sign1 =-1.;
        sign4 = 1.;
    }
    else if (almostEqual(chi[0]->evaluate(1.0), chi[3]->evaluate(1.0)))
    {
        sign1 =-1.;
        sign4 =-1.;
    }
    else if (almostEqual(chi[0]->evaluate(-1.0), chi[3]->evaluate(1.0)))
    {
        sign1 = 1.;
        sign4 =-1.;
    }
    else
    {
        std::println("TransfiniteQuadMap: chi1 and chi4 mismatch!");
        exit(EXIT_FAILURE);
    }
    if (almostEqual(chi[1]->evaluate(1.0), chi[2]->evaluate(1.0)))
    {
        sign2 = 1;
        sign3 = 1;
    }
    else if (almostEqual(chi[1]->evaluate(-1.0), chi[2]->evaluate(1.0)))
    {
        sign2 =-1;
        sign3 = 1;
    }
    else if (almostEqual(chi[1]->evaluate(-1.0), chi[2]->evaluate(-1.0)))
    {
        sign2 =-1;
        sign3 =-1;
    }
    else if (almostEqual(chi[1]->evaluate(1.0), chi[2]->evaluate(-1.0)))
    {
        sign2 = 1;
        sign3 =-1;
    }
    else
    {
        std::println("TransfiniteQuadMap: chi2 and chi3 mismatch!");
        exit(EXIT_FAILURE);
    }
    auto [x_1,  y_1] = chi[0]->evaluate(-sign1);
    auto [x_2,  y_2] = chi[0]->evaluate( sign1);
    auto [x_3,  y_3] = chi[2]->evaluate( sign3);
    auto [x_4,  y_4] = chi[2]->evaluate(-sign3);
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    arma::mat dxdx1(x1.size(), x2.size());
    arma::mat dydx1(x1.size(), x2.size());
    arma::mat dxdx2(x1.size(), x2.size());
    arma::mat dydx2(x1.size(), x2.size());
    for (size_t i = 0; i < x1.size(); i++)
    {
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_1, Ys_1] = chi[0]->derivative(sign1*x1i);
        auto [Xs_3, Ys_3] = chi[2]->derivative(sign3*x1i);
        for (size_t j = 0; j < x2.size(); j++)
        {
            double x2j = x2(j);
            auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
            auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
            auto [Xs_2, Ys_2] = chi[1]->derivative(sign2*x2j);
            auto [Xs_4, Ys_4] = chi[3]->derivative(sign4*x2j);
            dxdx1(i, j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
            dydx1(i, j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
            dxdx2(i, j) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
            dydx2(i, j) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
        }
    }
    return std::tie(dxdx1, dxdx2, dydx1, dydx2);
}