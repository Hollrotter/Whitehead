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
        b(i) =-arma::datum::tau*dot(Q, nC.col(i));
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
    dcp.fill(gamma_hat(0));
    for (size_t i = 0; i < nx; i++)
        dcp(i, 0) += 2*boost::math::chebyshev_clenshaw_recurrence(gamma_hat.memptr(), nx, xi(i));
    cL = 0;
    cM = 2./3;
    for (size_t k = 0; k < nx; k+=2)
        cL += 2./(1.-k*k)*gamma_hat(k);
    for (size_t k = 5; k < nx; k+=2)
        cM += 2./(4.-k*k)*gamma_hat(k);
}

void Airfoil::linearSolve()
{
    aerodynamicMatrix();
    arma::lu(L, U, P, A);
}

void Airfoil::linearEval()
{
    b(0) = 0;
    for (size_t i = 1; i < nx; i++)
        b(i) =-arma::datum::tau*alpha;
    gamma_hat = solve(trimatu(U), solve(trimatl(L), P*b));
}