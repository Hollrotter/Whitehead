#include "DVM.hpp"

/**
 * @brief 
 * 
 * @param _qdyn Dynamic pressure of the inflow.
 */
void DVM::dynamicPressure(double _qdyn)
{
    if (qdyn <= 0)
    {
        std::println("Dynamic pressure must be positive!");
        exit(EXIT_FAILURE);
    }
    qdyn = _qdyn;
}

void DVM::pitch(double _alpha)
{
    pitch(_alpha*arma::ones(1));
    con = 1;
}

void DVM::pitch(arma::vec _alpha)
{
    alpha = arma::datum::pi/180*_alpha;
    con   = alpha.size();
}

/**
 * @brief 
 * 
 */
void DVM::geometry()
{
    arma::vec x0 = x.head(nx);
    arma::vec x1 = x.tail(nx);
    arma::vec z0 = z.head(nx);
    arma::vec z1 = z.tail(nx);
    xg = (3*x0 + x1)/4;
    xC = (3*x1 + x0)/4;
    zg = (3*z0 + z1)/4;
    zC = (3*z1 + z0)/4;
}

/**
 * @brief 
 * 
 */
void DVM::dvm()
{
    switch(analysis)
    {
        case Analysis::linear:
            dvmSolve();
            dvmEval();
            break;
        case Analysis::nonlinear:
            dvmNonlinear();
            break;
    } 
}

/**
 * @brief 
 * 
 */
void DVM::dvmSolve()
{
    geometry();
    aerodynamicMatrix();
    arma::lu(L, U, P, A);
}

/**
 * @brief 
 * 
 */
void DVM::dvmEval()
{
    arma::vec w = diff(xC/c);
    arma::mat g = solve(trimatu(U), solve(trimatl(L), P*(repelem(w, 1, con) - repelem(alpha.t(), nx, 1))));
    postprocessing(g);
}

/**
 * @brief 
 * 
 */
void DVM::dvmNonlinear()
{
    z = c*camber(x/c);
    geometry();
    nC.zeros(2, nx);
    #pragma omp parallel for
    for (size_t i = 0; i < nx; i++)
        nC.col(i) = normalise(arma::vec::fixed<2>{-camber.diff(xC(i)/c), 1});
    aerodynamicMatrix();
    arma::mat Q = join_horiz(cos(alpha), sin(alpha));
    arma::mat g = solve(A,-(Q*nC).t());
    postprocessing(g);
}

/**
 * @brief 
 * 
 * @param filename File the data is written to.
 */
void DVM::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t i = 0; i < nx; i++)
        file << xg(i) << dcp.row(i);
    file.close();
}

/**
 * @brief 
 * 
 * @param g Circulation/inflow velocity calculated by Discrete-Vortex-Method.
 */
void DVM::postprocessing(arma::mat &g)
{
    cL = 2*sum(g).t()/c;
    switch(analysis)
    {
        case Analysis::linear:
            cM  = 2*g.t()*(c/4-xg)/pow(c, 2);
            dcp = 2*g/repelem(arma::diff(x), 1, con);
            break;
        case Analysis::nonlinear:
            cM  = 2*sum((c/4-xg*cos(alpha).t())%g).t()/pow(c, 2);
            dcp.zeros(nx, con);
            arma::vec alpha_panel =-atan(camber.diff(xC/c));
            for (size_t i = 0; i < nx; i++)
                dcp.row(i) = 2*cos(alpha + alpha_panel(i))%g.row(i)/sqrt(pow(x(i+1)-x(i), 2) + pow(z(i+1)-z(i), 2)); // Maybe wrong in Katz & Plotkin!
            break;
    }
}