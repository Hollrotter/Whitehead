#include "Wing.hpp"

Wing Wing::fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec xi_1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec xi_2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
        _h = Lagrange::covariantScaleFactors(xi_1, xi_2, _chi);
    return {_x, _y, _chi, _h};
}

Wing Wing::fromTransfiniteQuadMap(arma::mat _z, std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec xi_1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec xi_2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
        _h = Lagrange::covariantScaleFactors(xi_1, xi_2, _chi, _z);
    return {_x, _y, _z, _chi, _h};
}

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
            b.row(i+j*nx) =-4*arma::datum::pi*alpha;
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
    nonlinearSolve();
    nonlinearEval();
    postprocessing();
}

void Wing::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t i = 0; i < nx; i++, file << '\n')
        for (size_t j = 0; j < ny; j++, file << '\n')
        {
            file << x(i, j) << ' ' << y(i, j) << ' ' << z(i, j) << ' ' << mu(i, j);
            for (size_t k = 0; k < con; k++)
                file << ' ' << dcp(i, j, k);
        }
    file.close();
}

arma::mat Wing::calculateNormal()
{
    auto [dxdxi_1, dxdxi_2, dydxi_1, dydxi_2] = Lagrange::TransfiniteQuadMetrics(xi_1, xi_2, chi);
    arma::mat d2xdxi_12     = D1*dxdxi_1;
    arma::mat d2xdxi_1dxi_2 = dxdxi_1*D2.t();
    arma::mat d2xdxi_22     = dxdxi_2*D2.t();
    arma::mat d2ydxi_12     = D1*dydxi_1;
    arma::mat d2ydxi_1dxi_2 = D1*dydxi_2;
    arma::mat d2ydxi_22     = dydxi_2*D2.t();
    arma::mat dzdxi_1 = D1*zC;
    arma::mat dzdxi_2 = zC*D2.t();
    arma::mat d2zdxi_12 = D1*dzdxi_1;
    arma::mat d2zdxi_1dxi_2 = D1*dzdxi_2;
    arma::mat d2zdxi_22 = dzdxi_2*D2.t();

    arma::mat normal(nxy, 3);
    for (size_t j = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++)
        {
            double e_11 = e_c(i, j, 0);
            double e_12 = e_c(i, j, 1);
            double e_22 = e_c(i, j, 2);
            double e11  =  ec(i, j, 0);
            double e12  =  ec(i, j, 1);
            double e22  =  ec(i, j, 2);
            double e = e_11*e_22 - pow(e_12, 2);
            double sqrt_a = sqrt(e*(1 + e11*pow(dzdxi_1(i, j), 2) + 2*e12*dzdxi_1(i, j)*dzdxi_2(i, j) + e22*pow(dzdxi_2(i, j), 2)));
            normal.row(i+j*nx) = arma::rowvec::fixed<3>({dydxi_1(i, j)*dzdxi_2(i, j)-dzdxi_1(i, j)*dydxi_2(i, j),
                                                        dzdxi_1(i, j)*dxdxi_2(i, j)-dxdxi_1(i, j)*dzdxi_2(i, j),
                                                        dxdxi_1(i, j)*dydxi_2(i, j)-dydxi_1(i, j)*dxdxi_2(i, j)})/sqrt_a;
        }
    return normal;
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

void Wing::nonlinearSolve()
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

    arma::mat Q = join_horiz(cos(alpha), arma::zeros(con), sin(alpha)).t();
    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b.row(i+j*nx) =-2*arma::datum::tau*(nC.row(i+j*nx)*Q);

    arma::lu(L, U, P, A);
}

void Wing::nonlinearEval()
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
                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                {
                    double  t2 = boost::math::chebyshev_t(q, x2(j));
                    double dt2 = boost::math::chebyshev_t_prime(q, x2(j));
                    for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double  t1 = boost::math::chebyshev_t(p, x1(i));
                            double dt1 = boost::math::chebyshev_t_prime(p, x1(i));
                            mu(i, j)       += mu_hat(p+q*nx, 0)  * t1 * t2;
                            dcp.tube(i, j) += 2*mu_hat.row(p+q*nx) * (J11_inv(i, j)*dt1*t2 + J21_inv(i, j)*t1*dt2);
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
            arma::mat e = e_c.slice(0)%e_c.slice(2) - pow(e_c.slice(1), 2);
            arma::mat sqrt_a = sqrt(e%(1 + ec.slice(0)%pow(dzdx1, 2) + 2*ec.slice(1)%dzdx1%dzdx2 + ec.slice(2)%pow(dzdx2, 2)));
            for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                {
                    arma::vec::fixed<3> q_mu(arma::fill::zeros);
                    arma::vec::fixed<3> n = arma::vec::fixed<3>({J21(i, j)*dzdx2(i, j)-dzdx1(i, j)*J22(i, j),
                                                                 dzdx1(i, j)*J12(i, j)-J11(i, j)*dzdx2(i, j),
                                                                 J11(i, j)*J22(i, j)-J21(i, j)*J12(i, j)})/sqrt_a(i, j);
                    arma::mat::fixed<3, 2> J_red = {{J22(i, j)*n(2) - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - J21(i, j)*n(2)},
                                                    {n(0)*dzdx2(i, j) - J12(i, j)*n(2), J11(i, j)*n(2) - n(0)*dzdx1(i, j)},
                                                    {J12(i, j)*n(1) - n(0)*J22(i, j), n(0)*J21(i, j) - J11(i, j)*n(1)}};
                    J_red/=sqrt_a(i, j);
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
    auto [dxdx1_gl, dxdx2_gl, dydx1_gl, dydx2_gl] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);

    switch (analysis)
    {
        case Analysis::linear:
        {
            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    double detJ = dxdx1_gl(i, j)*dydx2_gl(i, j) - dxdx2_gl(i, j)*dydx1_gl(i, j);
                    double J11_inv = dydx2_gl(i, j)/detJ;
                    double J21_inv =-dydx1_gl(i, j)/detJ;
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
            arma::mat F(3, con);
            arma::mat M(3, con);
            arma::mat Q = join_horiz(cos(alpha), arma::zeros(con), sin(alpha));
            arma::mat Tx = Lagrange::interpolationMatrix(x1, x1_gl);
            arma::mat Ty = Lagrange::interpolationMatrix(x2, x2_gl);
            arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, z, x1_gl, x2_gl);
            arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
            arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
            arma::mat dzdx1_gl = D1_gl*z_gl;
            arma::mat dzdx2_gl = z_gl*D2_gl.t();
            arma::field<arma::mat> J_gl = {{dxdx1_gl, dxdx2_gl}, {dydx1_gl, dydx2_gl}};
            arma::cube e_c_gl = MetricCo(J_gl);
            arma::cube ec_gl  = MetricContra(e_c_gl);
            arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
            arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
            for (size_t j = 0; j < ny; j++)
                for (size_t i = 0; i < nx; i++)
                {
                    arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dydx1_gl(i, j)*dzdx2_gl(i, j)-dzdx1_gl(i, j)*dydx2_gl(i, j),
                                                                    dzdx1_gl(i, j)*dxdx2_gl(i, j)-dxdx1_gl(i, j)*dzdx2_gl(i, j),
                                                                    dxdx1_gl(i, j)*dydx2_gl(i, j)-dydx1_gl(i, j)*dxdx2_gl(i, j)})/sqrt_a(i, j);
                    arma::mat::fixed<3, 2> J_red = {{dydx2_gl(i, j)*n_gl(2) - n_gl(1)*dzdx2_gl(i, j),-(dydx1_gl(i, j)*n_gl(2) - n_gl(1)*dzdx1_gl(i, j))},
                                                    {-(dxdx2_gl(i, j)*n_gl(2) - n_gl(0)*dzdx2_gl(i, j)),dxdx1_gl(i, j)*n_gl(2) - n_gl(0)*dzdx1_gl(i, j)},
                                                    {dxdx2_gl(i, j)*n_gl(1) - n_gl(0)*dydx2_gl(i, j),-(dxdx1_gl(i, j)*n_gl(1) - n_gl(0)*dydx1_gl(i, j))}};
                    arma::vec::fixed<3> q_mu(arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++)
                    {
                        double T2 = boost::math::chebyshev_t(q, x2_gl(j));
                        double dT2dx2 = boost::math::chebyshev_t_prime(q, x2_gl(j));
                        for (size_t p = 0; p < nx; p++)
                        {
                            double T1 = boost::math::chebyshev_t(p, x1_gl(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, x1_gl(i));
                            arma::vec::fixed<2> dmudxi = {dT1dx1*T2, T1*dT2dx2};
                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                            q_mu += mu_hat(p+q*nx, 0) * gradmu;
                        }
                    }
                    arma::vec DCP = 2*Q*q_mu;
                    arma::vec r = {x_gl(i, j), y_gl(i, j), z_gl(i, j)};
                    area += gl_x[i].weight * gl_y[j].weight * sqrt_a(i, j);
                    F    += gl_x[i].weight * gl_y[j].weight * n_gl * DCP;
                    M    -= gl_x[i].weight * gl_y[j].weight * cross(n_gl * DCP, r);
                }
            lift   = F.row(2)*cos(alpha) - F.row(0)*sin(alpha);
            moment = M.row(1);
            break;
        }
    }
}