#include "Wing.hpp"

/**
 * @brief 
 * 
 * @param _qdyn Dynamic pressure of the inflow.
 */
void Wing::dynamicPressure(double _qdyn)
{
    qdyn = _qdyn;
    try
    {
        if (qdyn <= 0)
            throw std::runtime_error("Dynamic pressure must be positive!");
    }
    catch(const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief 
 * 
 * @param _alpha Pitch given as vector for multiple configurations.
 */
void Wing::pitch(arma::vec _alpha)
{
    alpha = arma::datum::pi/180*_alpha;
    con = alpha.size();
    lift.zeros(con);
    moment.zeros(con);
    dcp.zeros(nx, ny, con);
}

void Wing::checkMesh()
{
    std::println("Checking for negative volumes...");
    bool negativeVolumes = false;
    arma::mat detJ = J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0);
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
            if (detJ(i, j) < 0)
                negativeVolumes = true;
    if (negativeVolumes == false)
        std::println("No negative volumes found!");
    else
        std::println("Negative volumes were found!");
}

void Wing::linear()
{
    analysis = Analysis::linear;
    linearSolve();
    linearEval();
}

void Wing::nonlinear()
{
    analysis = Analysis::nonlinear;
    aerodynamicMatrix();
}

void Wing::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t m = 0; m < nx; m++, file << '\n')
        for (size_t n = 0; n < ny; n++, file << '\n')
        {
            file << x(m, n) << ' ' << y(m, n);
            for (size_t i = 0; i < con; i++)
                file << ' ' << dcp(m, n, i);
        }
    file.close();
}

void Wing::linearSolve()
{
    aerodynamicMatrix();
    
    for (size_t j = 1; j < ny-1; j++)
    {
        // BC west
        muBoundaryWest(j);
        // BC east
        muBoundaryEast(j);
    }
    for (size_t i = 1; i < nx-1; i++)
    {
        // BC south
        muBoundarySouth(i);
        // BC north
        muBoundaryNorth(i);
    }
    // BC south-west corner (i = 0, j = 0)
    if (chi[0]->curveType == CurveType::Interface && chi[3]->curveType == CurveType::Boundary)
        muBoundaryWest(0);
    else
        muBoundarySouth(0);

    // BC north-west corner (i = 0, j = ny-1)
    if (chi[2]->curveType == CurveType::Boundary  && chi[3]->curveType == CurveType::Interface)
        muBoundaryNorth(0);
    else
        muBoundaryWest(ny-1);

    // BC south-east corner (i = nx-1, j = 0)
    if (chi[0]->curveType == CurveType::Boundary  && chi[1]->curveType == CurveType::Interface)
        muBoundarySouth(nx-1);
    else
        muBoundaryEast(0);

    // BC north-east corner (i = nx-1, j = ny-1)
    if (chi[2]->curveType == CurveType::Interface && chi[1]->curveType == CurveType::Boundary)
        muBoundaryEast(ny-1);
    else
        muBoundaryNorth(nx-1);

    arma::lu(L, U, P, A);
}

void Wing::linearEval()
{
    mu_hat = solve(trimatu(U), solve(trimatl(L), P*b));
    postprocessing();
}

void Wing::postprocessing()
{
    arma::mat J11 = J(0, 0);
    arma::mat J12 = J(0, 1);
    arma::mat J21 = J(1, 0);
    arma::mat J22 = J(1, 1);
    arma::mat detJ = J11%J22 - J12%J21;
    arma::mat J11_inv = J22/detJ;
    arma::mat J12_inv =-J12/detJ;
    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
        for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
            {
                double T2 = boost::math::chebyshev_t(q, x2(j));
                double dT2dx2 = boost::math::chebyshev_t_prime(q, x2(j));
                for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                {
                    double T1 = boost::math::chebyshev_t(p, x1(i));
                    double dT1dx1 = boost::math::chebyshev_t_prime(p, x1(i));
                    mu(i, j)       += mu_hat(p+q*nx, 0)  * T1 * T2;
                    dcp.tube(i, j) += mu_hat.row(p+q*nx) * (J11_inv(i, j)*dT1dx1*T2 + J12_inv(i, j)*T1*dT2dx2);
                }
            }
}