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
    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b.row(i+j*nx) =-4*arma::datum::pi*alpha; // Must be corrected later!
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
    postprocessing();
}

void Wing::nonlinear()
{
    analysis = Analysis::nonlinear;
    nC.zeros(nxy, 3);
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

    arma::mat Q = join_horiz(cos(alpha), arma::zeros(con), sin(alpha)).t();
    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b.row(i+j*nx) =-2*arma::datum::tau*(nC.row(i+j*nx)*Q);
    mu_hat = solve(A, b);
    postprocessing();
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
}

void Wing::postprocessing()
{
    mu.zeros();
    dcp.zeros();
    arma::mat J11 = J(0, 0);
    arma::mat J12 = J(0, 1);
    arma::mat J21 = J(1, 0);
    arma::mat J22 = J(1, 1);
    switch (analysis)
    {
        case Analysis::linear:
        {
            arma::mat detJ = J11%J22 - J12%J21;
            arma::mat J11_inv = J22/detJ;
            arma::mat J21_inv =-J21/detJ;
            for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double T2 = boost::math::chebyshev_t(q, x2(j));
                        double dT2dx2 = boost::math::chebyshev_t_prime(q, x2(j));
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1 = boost::math::chebyshev_t(p, x1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, x1(i));
                            mu(i, j)       += mu_hat(p+q*nx, 0)  * T1 * T2;
                            dcp.tube(i, j) += 2*mu_hat.row(p+q*nx) * (J11_inv(i, j)*dT1dx1*T2 + J21_inv(i, j)*T1*dT2dx2);
                        }
                    }
            break;
        }
        case Analysis::nonlinear:
        {
            arma::mat d1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
            arma::mat d2 = Chebyshev::derivativeMatrix(x2, Derivative::first);
            arma::mat dzdx1 = d1 * z;
            arma::mat dzdx2 = z * d2.t();
            arma::mat Q = join_horiz(cos(alpha), arma::zeros(con), sin(alpha));
            for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                {
                    arma::vec::fixed<3> q_mu(arma::fill::zeros);
                    arma::vec::fixed<3> n = {J21(i, j)*dzdx2(i, j)-dzdx1(i, j)*J22(i, j),
                                             dzdx1(i, j)*J12(i, j)-J11(i, j)*dzdx2(i, j),
                                             J11(i, j)*J22(i, j)-J21(i, j)*J12(i, j)};
                    arma::mat::fixed<3, 2> J_red = {{J22(i, j)*n(2) - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - J21(i, j)*n(2)},
                                                    {n(0)*dzdx2(i, j) - J12(i, j)*n(2), J11(i, j)*n(2) - n(0)*dzdx1(i, j)},
                                                    {J12(i, j)*n(1) - n(0)*J22(i, j), n(0)*J21(i, j) - J11(i, j)*n(1)}};
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double T2 = boost::math::chebyshev_t(q, x2(j));
                        double dT2dx2 = boost::math::chebyshev_t_prime(q, x2(j));
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1 = boost::math::chebyshev_t(p, x1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, x1(i));
                            arma::vec::fixed<2> dmudxi = {dT1dx1*T2, T1*dT2dx2};
                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                            mu(i, j) += mu_hat(p+q*nx, 0) * T1 * T2;
                            q_mu     += mu_hat(p+q*nx, 0) * gradmu;
                        }
                    }
                    dcp.tube(i, j) = 2*Q*q_mu;
                }
            break;
        }
    }
    area = 0;
    lift.zeros();
    moment.zeros();
    std::vector<fastgl::QuadPair> gl_x(nx), gl_y(ny);
    arma::vec x1_gl(nx), x2_gl(ny);
    for (size_t i = 0; i < nx; i++)
    {
        gl_x[i] = fastgl::GLPair(nx, i+1);
        x1_gl(i) =-gl_x[i].x();
    }
    for (size_t j = 0; j < ny; j++)
    {
        gl_y[j] = fastgl::GLPair(ny, j+1);
        x2_gl(j) =-gl_y[j].x();
    }
    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);

    switch (analysis)
    {
        case Analysis::linear:
        {
            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    double detJ = dx_gldx1(i, j)*dy_gldx2(i, j) - dx_gldx2(i, j)*dy_gldx1(i, j);
                    double J11_inv = dy_gldx2(i, j)/detJ;
                    double J21_inv =-dy_gldx1(i, j)/detJ;
                    arma::vec DCP(con, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++)
                    {
                        double T2 = boost::math::chebyshev_t(q, x2_gl(j));
                        double dT2dx2 = boost::math::chebyshev_t_prime(q, x2_gl(j));
                        for (size_t p = 0; p < nx; p++)
                        {
                            double T1 = boost::math::chebyshev_t(p, x1_gl(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, x1_gl(i));
                            DCP += 2*mu_hat.row(p+q*nx) * (J11_inv*dT1dx1*T2 + J21_inv*T1*dT2dx2);
                        }
                    }
                    area   += gl_x[i].weight * gl_y[j].weight * detJ;
                    lift   += gl_x[i].weight * gl_y[j].weight * DCP * detJ;
                    moment -= gl_x[i].weight * gl_y[j].weight * DCP * detJ * x_gl(i, j);
                }
            break;
        }
        case Analysis::nonlinear:
        {
            break;
        }
    }
}