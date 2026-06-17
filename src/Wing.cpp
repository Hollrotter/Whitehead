#include "Wing.hpp"

Wing Wing::fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec _x1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec _x2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
        _h = Lagrange::covariantScaleFactors(_x1, _x2, _chi);
    return {_chi, _x, _y, _h};
}

Wing Wing::fromTransfiniteQuadMap(arma::mat _z, std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec _x1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec _x2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
        _h = Lagrange::covariantScaleFactors(_x1, _x2, _chi, _z);
    return {_chi, _x, _y, _z, _h};
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
            file << x(i, j) << ' ' << y(i, j) << ' ' << z(i, j) << ' ' << mu(i, j) << ' ' << dcp(i, j);
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

    arma::mat normal(nxy, 3, arma::fill::none);
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

    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b(i+j*nx) =-2*arma::datum::tau*alpha;

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

    arma::vec Q = {cos(alpha), 0, sin(alpha)};
    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b(i+j*nx) =-2*arma::datum::tau*dot(nC.row(i+j*nx), Q);

    arma::lu(L, U, P, A);
}

void Wing::nonlinearEval()
{
    mu_hat = solve(trimatu(U), solve(trimatl(L), P*b));
}

void Wing::postprocessing()
{
    arma::mat J11 = J(0, 0);
    arma::mat J12 = J(0, 1);
    arma::mat J21 = J(1, 0);
    arma::mat J22 = J(1, 1);
    arma::mat MU_0 = reshape(mu_hat, nx, ny);
    auto [MU_1, MU_2] = Chebyshev::DerivativeCoefficients(MU_0);
    arma::vec MU_0_y(ny, arma::fill::none);
    arma::vec MU_1_y(ny, arma::fill::none);
    arma::vec MU_2_y(ny, arma::fill::none);
    switch (analysis)
    {
        case Analysis::linear:
        {
            arma::mat detJ = J11%J22 - J12%J21;
            arma::mat J11_inv = J22/detJ;
            arma::mat J21_inv =-J21/detJ;
            for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
            {
                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                {
                    MU_0_y(q) = MU_0(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_0.colptr(q), nx, x1(i));
                    MU_1_y(q) = MU_1(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1.colptr(q), nx, x1(i));
                    MU_2_y(q) = MU_2(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2.colptr(q), nx, x1(i));
                }
                for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                {
                    mu(i, j)  = MU_0_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_0_y.memptr(), ny, x2(j));
                    dcp(i, j) = 2*(J11_inv(i, j)*(MU_1_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1_y.memptr(), ny, x2(j)))
                                 + J21_inv(i, j)*(MU_2_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2_y.memptr(), ny, x2(j))));
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
            arma::vec Q = {cos(alpha), 0, sin(alpha)};
            arma::mat e = e_c.slice(0)%e_c.slice(2) - pow(e_c.slice(1), 2);
            arma::mat sqrt_a = sqrt(e%(1 + ec.slice(0)%pow(dzdx1, 2) + 2*ec.slice(1)%dzdx1%dzdx2 + ec.slice(2)%pow(dzdx2, 2)));
            for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
            {
                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                {
                    MU_0_y(q) = MU_0(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_0.colptr(q), nx, x1(i));
                    MU_1_y(q) = MU_1(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1.colptr(q), nx, x1(i));
                    MU_2_y(q) = MU_2(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2.colptr(q), nx, x1(i));
                }
                for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                {
                    arma::vec::fixed<3> n = arma::vec::fixed<3>({J21(i, j)*dzdx2(i, j)-dzdx1(i, j)*J22(i, j),
                                                                 dzdx1(i, j)*J12(i, j)-J11(i, j)*dzdx2(i, j),
                                                                 J11(i, j)*J22(i, j)-J21(i, j)*J12(i, j)})/sqrt_a(i, j);
                    arma::mat::fixed<3, 2> J_red = {{J22(i, j)*n(2) - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - J21(i, j)*n(2)},
                                                    {n(0)*dzdx2(i, j) - J12(i, j)*n(2), J11(i, j)*n(2) - n(0)*dzdx1(i, j)},
                                                    {J12(i, j)*n(1) - n(0)*J22(i, j), n(0)*J21(i, j) - J11(i, j)*n(1)}};
                    J_red/=sqrt_a(i, j);
                    mu(i, j) = MU_0_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_0_y.memptr(), ny, x2(j));
                    arma::vec::fixed<2> dmudxi = {MU_1_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1_y.memptr(), ny, x2(j)),
                                                  MU_2_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2_y.memptr(), ny, x2(j))};
                    arma::vec::fixed<3> q_mu = J_red*dmudxi;
                    dcp(i, j) = 2*dot(Q, q_mu);
                }
            }
            break;
        }
        default:
            std::println("Only linear and nonlinear analysis are implemented for Wing!");
            exit(EXIT_FAILURE);
    }
    area   = 0;
    lift   = 0;
    moment = 0;
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
            for (size_t i = 0; i < nx; i++)
            {
                for (size_t q = 0; q < ny; q++)
                {
                    MU_1_y(q) = MU_1(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1.colptr(q), nx, x1_gl(i));
                    MU_2_y(q) = MU_2(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2.colptr(q), nx, x1_gl(i));
                }
                for (size_t j = 0; j < ny; j++)
                {
                    double detJ = dxdx1_gl(i, j)*dydx2_gl(i, j) - dxdx2_gl(i, j)*dydx1_gl(i, j);
                    double J11_inv = dydx2_gl(i, j)/detJ;
                    double J21_inv =-dydx1_gl(i, j)/detJ;
                    double DCP = 2*(J11_inv*(MU_1_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1_y.memptr(), ny, x2_gl(j)))
                                  + J21_inv*(MU_2_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2_y.memptr(), ny, x2_gl(j))));
                    area   += gl_x[i].weight * gl_y[j].weight * detJ;
                    lift   += gl_x[i].weight * gl_y[j].weight * DCP * detJ;
                    moment -= gl_x[i].weight * gl_y[j].weight * DCP * detJ * x_gl(i, j);
                }
            }
            break;
        }
        case Analysis::nonlinear:
        {
            arma::vec F(3);
            arma::vec M(3);
            arma::vec Q = {cos(alpha), 0, sin(alpha)};
            arma::mat Tx_gl = Lagrange::interpolationMatrix(x1, x1_gl);
            arma::mat Ty_gl = Lagrange::interpolationMatrix(x2, x2_gl);
            arma::mat z_gl  = Lagrange::interpolation2D(Tx_gl, Ty_gl, z, x1_gl, x2_gl);
            arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
            arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
            arma::mat dzdx1_gl = D1_gl*z_gl;
            arma::mat dzdx2_gl = z_gl*D2_gl.t();
            arma::field<arma::mat> J_gl = {{dxdx1_gl, dxdx2_gl}, {dydx1_gl, dydx2_gl}};
            arma::cube e_c_gl = MetricCo(J_gl);
            arma::cube ec_gl  = MetricContra(e_c_gl);
            arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
            arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
            for (size_t i = 0; i < nx; i++)
            {
                for (size_t q = 0; q < ny; q++)
                {
                    MU_1_y(q) = MU_1(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1.colptr(q), nx, x1_gl(i));
                    MU_2_y(q) = MU_2(0, q)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2.colptr(q), nx, x1_gl(i));
                }
                for (size_t j = 0; j < ny; j++)
                {
                    arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dydx1_gl(i, j)*dzdx2_gl(i, j)-dzdx1_gl(i, j)*dydx2_gl(i, j),
                                                                    dzdx1_gl(i, j)*dxdx2_gl(i, j)-dxdx1_gl(i, j)*dzdx2_gl(i, j),
                                                                    dxdx1_gl(i, j)*dydx2_gl(i, j)-dydx1_gl(i, j)*dxdx2_gl(i, j)})/sqrt_a(i, j);
                    arma::mat::fixed<3, 2> J_red = {{dydx2_gl(i, j)*n_gl(2)-n_gl(1)*dzdx2_gl(i, j), n_gl(1)*dzdx1_gl(i, j)-dydx1_gl(i, j)*n_gl(2)},
                                                    {dzdx2_gl(i, j)*n_gl(0)-n_gl(2)*dxdx2_gl(i, j), n_gl(2)*dxdx1_gl(i, j)-dzdx1_gl(i, j)*n_gl(0)},
                                                    {dxdx2_gl(i, j)*n_gl(1)-n_gl(0)*dydx2_gl(i, j), n_gl(0)*dydx1_gl(i, j)-dxdx1_gl(i, j)*n_gl(1)}};

                    arma::vec::fixed<2> dmudxi = {MU_1_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_1_y.memptr(), ny, x2_gl(j)),
                                                  MU_2_y(0)/2 + boost::math::chebyshev_clenshaw_recurrence(MU_2_y.memptr(), ny, x2_gl(j))};
                    arma::vec::fixed<3> q_mu = J_red*dmudxi;
                    double DCP = 2*dot(Q, q_mu);
                    arma::vec r = {x_gl(i, j), y_gl(i, j), z_gl(i, j)};
                    area += gl_x[i].weight * gl_y[j].weight * sqrt_a(i, j);
                    F    += gl_x[i].weight * gl_y[j].weight * n_gl * DCP;
                    M    -= gl_x[i].weight * gl_y[j].weight * cross(n_gl * DCP, r);
                }
            }
            lift   = F(2)*cos(alpha) - F(0)*sin(alpha);
            moment = M(1);
            break;
        }
        default:
            std::println("Only linear and nonlinear analysis are implemented for Wing!");
            exit(EXIT_FAILURE);
    }
}