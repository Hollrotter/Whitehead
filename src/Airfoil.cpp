#include "Airfoil.hpp"

Airfoil Airfoil::fromLagrangeCurveInterpolant(Lagrange::CurveInterpolant* _chi)
{
    arma::vec _xi = Chebyshev::gauss(_chi->getNodes().size());
    auto [_x, _z] = _chi->evaluate(_xi);
    return {_x, _z, _chi};
}

/**
 * @brief 
 * 
 * @param _qdyn Dynamic pressure of the inflow.
 */
void Airfoil::dynamicPressure(double _qdyn)
{
    if (qdyn <= 0)
    {
        std::println("Dynamic pressure must be positive!");
        exit(EXIT_FAILURE);
    }
    qdyn = _qdyn;
}

void Airfoil::pitch(double _alpha)
{
    alpha = arma::datum::pi/180*_alpha;
}

void Airfoil::linear()
{
    analysis = Analysis::linear;
    linearSolve();
    linearEval();
    postprocessing();
}

void Airfoil::nonlinear()
{
    analysis = Analysis::nonlinear;
    aerodynamicMatrix();
    arma::vec Q = {cos(alpha), sin(alpha)};
    b(0) = 0;
    for (size_t i = 1; i < nx; i++)
        b(i) = 2*dot(Q, nC.col(i));
    gamma_hat = solve(A, b);
    postprocessing();
}

/**
 * @brief 
 * 
 * @param filename File the data is written to.
 */
void Airfoil::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t i = 0; i < nx; i++)
        file << x(i) << ' ' << dcp(i) << '\n';
    file.close();
}

void Airfoil::postprocessing()
{
    dcp.zeros(nx);
    if (analysis == Analysis::linear)
        for (size_t i = 0; i < nx; i++)
        {
            double w = 1; // W_j(xi): Chebyshev Polynomial of the fourth kind
            double wp1 = 2*xi(i) + 1; // W_j+1(xi): Next Chebyshev Polynomial
            for (size_t j = 0; j < nx; j++)
            {
                dcp(i) += gamma_hat(j)*w;
                std::swap(w, wp1);
                wp1 = boost::math::chebyshev_next(xi(i), w, wp1);
            }
            dcp(i) *= 2*sqrt((1-xi(i))/(1+xi(i)));
        }
    else
    {
        auto [dx, dz] = chi->derivative(xi);
        arma::vec alpha_panel =-atan(dz/dx);
        for (size_t i = 0; i < nx; i++)
        {
            double w = 1; // W_j(xi): Chebyshev Polynomial of the fourth kind
            double wp1 = 2*xi(i) + 1; // W_j+1(xi): Next Chebyshev Polynomial
            for (size_t j = 0; j < nx; j++)
            {
                dcp(i) += gamma_hat(j)*w;
                std::swap(w, wp1);
                wp1 = boost::math::chebyshev_next(xi(i), w, wp1);
            }
            dcp(i) *= 2*cos(alpha + alpha_panel(i))*sqrt((1-xi(i))/(1+xi(i)));
        }
    }
    cL = arma::datum::pi * gamma_hat(0);
    cM = 0;
    for (size_t k = 0; k < nx; k++)
    {
        double x_k = cos(arma::datum::pi*(nx-1)/(nx+0.5));
        double W = 1;
        double Wp1 = 2*x_k + 1;
        for (size_t j = 0; j < nx; j++)
        {
            cM += arma::datum::pi*(1-x_k)/(nx+0.5) * x_k*gamma_hat(j)*W;
            std::swap(W, Wp1);
            Wp1 = boost::math::chebyshev_next(x_k, W, Wp1);
        }
    }
}

void Airfoil::linearSolve()
{
    aerodynamicMatrix();
    arma::lu(L, U, P, A);
}

void Airfoil::linearEval()
{
    for (size_t i = 0; i < nx; i++)
        b(i) = 2*(alpha - camber.diff(x(i)/c));
    gamma_hat = solve(trimatu(U), solve(trimatl(L), P*b));
}