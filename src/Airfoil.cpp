#include "Airfoil.hpp"

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
    pitch(_alpha*arma::ones(1));
    con = 1;
}

void Airfoil::pitch(arma::vec _alpha)
{
    alpha = arma::datum::pi/180*_alpha;
    con   = alpha.size();
}

void Airfoil::linear()
{
    analysis = Analysis::linear;
    linearSolve();
    linearEval();
}

void Airfoil::nonlinear()
{
    analysis = Analysis::nonlinear;
    aerodynamicMatrix();
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
        file << x(i) << dcp.row(i);
    file.close();
}

void Airfoil::postprocessing()
{
    dcp.zeros(nx, con);
    for (size_t k = 0; k < nx; k++)
        dcp += 2*Chebyshev::Polynomial(k, xi)*gamma_hat.row(k);
    cL.zeros(con);
    cM.zeros(con);
    cM.fill(2./3);
    for (size_t k = 0; k < nx; k+=2)
        cL += 2./(1.-k*k)*gamma_hat.row(k);
    for (size_t k = 5; k < nx; k+=2)
        cM += 2./(4.-k*k)*gamma_hat.row(k);
}

void Airfoil::linearSolve()
{
    aerodynamicMatrix();
    arma::lu(L, U, P, A);
}

void Airfoil::linearEval()
{
    gamma_hat = solve(trimatu(U), solve(trimatl(L), P*b));
    postprocessing();
}