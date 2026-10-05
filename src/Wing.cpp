#include "Wing.hpp"

Wing Wing::fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec _x1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec _x2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    return {_chi, _x, _y};
}

Wing Wing::fromTransfiniteQuadMap(arma::mat _z, std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
    arma::vec _x1 = Chebyshev::gauss(_chi[0]->getNodes().size());
    arma::vec _x2 = Chebyshev::gauss(_chi[1]->getNodes().size());
    return {_chi, _x, _y, _z};
}

void Wing::checkMesh() const
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
    if (mu.westBC == BC::None && mu.eastBC == BC::None && chi[3]->curveType == CurveType::Boundary && chi[1]->curveType == CurveType::Boundary)
        phi1.reset(new ChebyshevU(nx, mx));
    else if (mu.westBC == BC::None && chi[3]->curveType == CurveType::Boundary)
        phi1.reset(new JacobiBeta(nx, mx));
    else if (mu.eastBC == BC::None && chi[1]->curveType == CurveType::Boundary)
        phi1.reset(new JacobiAlpha(nx, mx));
    else
        phi1.reset(new ChebyshevT(nx, mx));
    if (mu.southBC == BC::None && mu.northBC == BC::None && chi[0]->curveType == CurveType::Boundary && chi[2]->curveType == CurveType::Boundary)
        phi2.reset(new ChebyshevU(ny, my));
    else if (mu.southBC == BC::None && chi[0]->curveType == CurveType::Boundary)
        phi2.reset(new JacobiBeta(ny, my));
    else if (mu.northBC == BC::None && chi[2]->curveType == CurveType::Boundary)
        phi2.reset(new JacobiAlpha(ny, my));
    else
        phi2.reset(new ChebyshevT(ny, my));
    init();
    linearSolve();
    linearEval();
    postprocessing();
}

void Wing::nonlinear()
{
    analysis = Analysis::nonlinear;
    if (mu.westBC == BC::None && mu.eastBC == BC::None && chi[3]->curveType == CurveType::Boundary && chi[1]->curveType == CurveType::Boundary)
        phi1.reset(new ChebyshevU(nx, mx));
    else if (mu.westBC == BC::None && chi[3]->curveType == CurveType::Boundary)
        phi1.reset(new JacobiBeta(nx, mx));
    else if (mu.eastBC == BC::None && chi[1]->curveType == CurveType::Boundary)
        phi1.reset(new JacobiAlpha(nx, mx));
    else
        phi1.reset(new ChebyshevT(nx, mx));
    if (mu.southBC == BC::None && mu.northBC == BC::None && chi[0]->curveType == CurveType::Boundary && chi[2]->curveType == CurveType::Boundary)
        phi2.reset(new ChebyshevU(ny, my));
    else if (mu.southBC == BC::None && chi[0]->curveType == CurveType::Boundary)
        phi2.reset(new JacobiBeta(ny, my));
    else if (mu.northBC == BC::None && chi[2]->curveType == CurveType::Boundary)
        phi2.reset(new JacobiAlpha(ny, my));
    else
        phi2.reset(new ChebyshevT(ny, my));
    init();
    nonlinearSolve();
    nonlinearEval();
    postprocessing();
}

double Wing::get_area()
{
    if (areaComputed == false)
    {
        area = 0;
        arma::mat J11 = J(0, 0);
        arma::mat J12 = J(0, 1);
        arma::mat J21 = J(1, 0);
        arma::mat J22 = J(1, 1);
        std::vector<fastgl::QuadPair> gl_x(nx), gl_y(ny);
        arma::vec x1_gl(nx, arma::fill::none), x2_gl(ny, arma::fill::none);
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
                #pragma omp parallel for reduction(+:area)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                        area += gl_x[i].weight * gl_y[j].weight * (dxdx1_gl(i, j)*dydx2_gl(i, j) - dxdx2_gl(i, j)*dydx1_gl(i, j));
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec F(3, arma::fill::zeros), M(3, arma::fill::zeros);
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
                #pragma omp parallel for reduction(+:area)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                        area += gl_x[i].weight * gl_y[j].weight * sqrt_a(i, j);
                break;
            }
            default:
                std::println("Only linear and nonlinear analysis are implemented for Wing!");
                exit(EXIT_FAILURE);
        }
        areaComputed = true;
    }
    return area;
}

void Wing::output(std::string filename) const
{
    std::ofstream file(filename);
    for (size_t i = 0; i < nx; i++, file << '\n')
        for (size_t j = 0; j < ny; j++, file << '\n')
            file << x(i, j) << ' ' << y(i, j) << ' ' << z(i, j) << ' ' << mu(i, j) << ' ' << dcp(i, j);
    file.close();
}

void Wing::init()
{
    A.zeros();
    b.zeros();
    xi_1 = phi1->xi;
    xi_2 = phi2->xi;
    PHI1.col(0)  = phi1->constant(nx);
    PHI2.col(0)  = phi2->constant(ny);
    dPHI1.col(0) = phi1->constantDerivative(nx);
    dPHI2.col(0) = phi2->constantDerivative(ny);
    PHI1.col(1)  = phi1->linear(xi_1);
    PHI2.col(1)  = phi2->linear(xi_2);
    dPHI1.col(1) = phi1->linearDerivative(nx);
    dPHI2.col(1) = phi2->linearDerivative(ny);
    for (size_t p = 1; p < nx-1; p++)
    {
        phi1->next(p, xi_1, PHI1);
        phi1->nextDerivative(p, xi_1, PHI1, dPHI1);
    }
    for (size_t q = 1; q < ny-1; q++)
    {
        phi2->next(q, xi_2, PHI2);
        phi2->nextDerivative(q, xi_2, PHI2, dPHI2);
    }
    D1 = Lagrange::derivativeMatrix(xi_1);
    D2 = Lagrange::derivativeMatrix(xi_2);
    std::tie(xC, yC) = Lagrange::TransfiniteQuadMap(xi_1, xi_2, chi);
    if (analysis == Analysis::linear)
        std::tie(h_2s2_south, h_2s1_south, h_1s1_east, h_1s2_east, h_2s2_north, h_2s1_north, h_1s1_west, h_1s2_west)
            = Lagrange::covariantScaleFactors(xi_1, xi_2, chi);
    else
    {
        Tx = Lagrange::interpolationMatrix(x1, xi_1);
        Ty = Lagrange::interpolationMatrix(x2, xi_2);
        zC = Lagrange::interpolation2D(Tx, Ty, z, xi_1, xi_2);
        std::tie(h_2s2_south, h_2s1_south, h_1s1_east, h_1s2_east, h_2s2_north, h_2s1_north, h_1s1_west, h_1s2_west)
            = Lagrange::covariantScaleFactors(xi_1, xi_2, chi, z, D1, D2);
        nC = calculateNormal();
    }
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
            double e11  =  ec(i, j, 0);
            double e12  =  ec(i, j, 1);
            double e22  =  ec(i, j, 2);
            double e = e_c(i, j, 0)*e_c(i, j, 2) - pow(e_c(i, j, 1), 2);
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
    if (phi1->basis == Basis::T || phi1->basis == Basis::PA)
        muBoundaryWest(0);
    else
        muBoundarySouth(0);

    // BC north-west corner (i = 0, j = ny-1)
    if (phi2->basis == Basis::T || phi2->basis == Basis::PB)
        muBoundaryNorth(0);
    else
        muBoundaryWest(ny-1);

    // BC south-east corner (i = nx-1, j = 0)
    if (phi2->basis == Basis::T || phi2->basis == Basis::PA)
        muBoundarySouth(nx-1);
    else
        muBoundaryEast(0);

    // BC north-east corner (i = nx-1, j = ny-1)
    if (phi1->basis == Basis::T || phi1->basis == Basis::PB)
        muBoundaryEast(ny-1);
    else
        muBoundaryNorth(nx-1);

    size_t i_min = 0, i_max = nx, j_min = 0, j_max = ny;
    if (phi1->basis == Basis::T || phi1->basis == Basis::PA)
        i_min = 1;
    if (phi1->basis == Basis::T || phi1->basis == Basis::PB)
        i_max = nx-1;
    if (phi2->basis == Basis::T || phi2->basis == Basis::PA)
        j_min = 1;
    if (phi2->basis == Basis::T || phi2->basis == Basis::PB)
        j_max = ny-1;
    for (size_t j = j_min; j < j_max; j++)
        for (size_t i = i_min; i < i_max; i++)
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
    if (phi1->basis == Basis::T || phi1->basis == Basis::PA)
        muBoundaryWest(0);
    else
        muBoundarySouth(0);

    // BC north-west corner (i = 0, j = ny-1)
    if (phi2->basis == Basis::T || phi2->basis == Basis::PB)
        muBoundaryNorth(0);
    else
        muBoundaryWest(ny-1);

    // BC south-east corner (i = nx-1, j = 0)
    if (phi2->basis == Basis::T || phi2->basis == Basis::PA)
        muBoundarySouth(nx-1);
    else
        muBoundaryEast(0);

    // BC north-east corner (i = nx-1, j = ny-1)
    if (phi1->basis == Basis::T || phi1->basis == Basis::PB)
        muBoundaryEast(ny-1);
    else
        muBoundaryNorth(nx-1);

    arma::vec Q = {cos(alpha), 0, sin(alpha)};
    size_t i_min = 0, i_max = nx, j_min = 0, j_max = ny;
    if (phi1->basis == Basis::T || phi1->basis == Basis::PA)
        i_min = 1;
    if (phi1->basis == Basis::T || phi1->basis == Basis::PB)
        i_max = nx-1;
    if (phi2->basis == Basis::T || phi2->basis == Basis::PA)
        j_min = 1;
    if (phi2->basis == Basis::T || phi2->basis == Basis::PB)
        j_max = ny-1;
    for (size_t j = j_min; j < j_max; j++)
        for (size_t i = i_min; i < i_max; i++)
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
    if (phi1->basis == Basis::T && phi2->basis == Basis::T)
    {
        arma::mat MU_0 = reshape(mu_hat, nx, ny);
        auto [MU_1, MU_2] = Chebyshev::DerivativeCoefficients(MU_0);
        switch (analysis)
        {
            case Analysis::linear:
            {
                arma::mat detJ = J11%J22 - J12%J21;
                arma::mat J11_inv = J22/detJ;
                arma::mat J21_inv =-J21/detJ;
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                {
                    arma::vec MU_0_y(ny, arma::fill::none), MU_1_y(ny, arma::fill::none), MU_2_y(ny, arma::fill::none);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomials in 2-direction
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
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                {
                    arma::vec MU_0_y(ny, arma::fill::none), MU_1_y(ny, arma::fill::none), MU_2_y(ny, arma::fill::none);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomials in 2-direction
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
                        arma::mat::fixed<3, 2> J_red = {{n(2)*J22(i, j)   - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - n(2)*J21(i, j)},
                                                        {n(0)*dzdx2(i, j) - n(2)*J12(i, j),   n(2)*J11(i, j)   - n(0)*dzdx1(i, j)},
                                                        {n(1)*J12(i, j)   - n(0)*J22(i, j),   n(0)*J21(i, j)   - n(1)*J11(i, j)}};
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
        arma::vec x1_gl = phi1->xg;
        arma::vec x2_gl = phi2->xg;
        arma::vec w1_gl = phi1->wg;
        arma::vec w2_gl = phi2->wg;
        auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
        auto [dxdx1_gl, dxdx2_gl, dydx1_gl, dydx2_gl] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);

        switch (analysis)
        {
            case Analysis::linear:
            {
                #pragma omp parallel for reduction(+:area) reduction(+:lift) reduction(+:moment)
                for (size_t i = 0; i < nx; i++)
                {
                    arma::vec MU_1_y(ny, arma::fill::none);
                    arma::vec MU_2_y(ny, arma::fill::none);
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
                        double dA = w1_gl(i) * w2_gl(j) * detJ;
                        area   += dA;
                        lift   += dA * DCP;
                        moment -= dA * DCP * x_gl(i, j);
                    }
                }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec F(3, arma::fill::zeros), M(3, arma::fill::zeros);
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
                #pragma omp parallel for reduction(+:area) reduction(+:F) reduction(+:M)
                for (size_t i = 0; i < nx; i++)
                {
                    arma::vec MU_1_y(ny, arma::fill::none);
                    arma::vec MU_2_y(ny, arma::fill::none);
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
                        area += w1_gl(i) * w2_gl(j) * sqrt_a(i, j);
                        F    += w1_gl(i) * w2_gl(j) * n_gl * DCP;
                        M    -= w1_gl(i) * w2_gl(j) * cross(n_gl * DCP, r);
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
        areaComputed = true;
    }
    else if (phi1->basis == Basis::T)
    {
        mu.zeros();
        dcp.zeros();
        switch (analysis)
        {
            case Analysis::linear:
            {
                arma::mat detJ = J11%J22 - J12%J21;
                arma::mat J11_inv = J22/detJ;
                arma::mat J21_inv =-J21/detJ;
                arma::vec  w2 = phi2->weightFunction(x2);
                arma::vec dw2 = phi2->weightFunctionDerivative(x2);
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                double  psi2 = w2(j)*t2;
                                double dpsi2 = w2(j)*dt2 + dw2(j)*t2;
                                mu(i, j)  +=   mu_hat(p+q*nx) * t1*psi2;
                                dcp(i, j) += 2*mu_hat(p+q*nx) * (J11_inv(i, j)*dt1*psi2 + J21_inv(i, j)*t1*dpsi2);
                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec  w2 = phi2->weightFunction(x2);
                arma::vec dw2 = phi2->weightFunctionDerivative(x2);
                arma::mat d1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
                arma::mat d2 = Chebyshev::derivativeMatrix(x2, Derivative::first);
                arma::mat dzdx1 = d1 * z;
                arma::mat dzdx2 = z * d2.t();
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat e = e_c.slice(0)%e_c.slice(2) - pow(e_c.slice(1), 2);
                arma::mat sqrt_a = sqrt(e%(1 + ec.slice(0)%pow(dzdx1, 2) + 2*ec.slice(1)%dzdx1%dzdx2 + ec.slice(2)%pow(dzdx2, 2)));
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        arma::vec::fixed<3> n = arma::vec::fixed<3>({J21(i, j)*dzdx2(i, j)-J22(i, j)*dzdx1(i, j),
                                                                     J12(i, j)*dzdx1(i, j)-J11(i, j)*dzdx2(i, j),
                                                                     J11(i, j)*J22(i, j)  -J21(i, j)*J12(i, j)})/sqrt_a(i, j);
                        arma::mat::fixed<3, 2> J_red = {{n(2)*J22(i, j)   - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - n(2)*J21(i, j)},
                                                        {n(0)*dzdx2(i, j) - n(2)*J12(i, j),   n(2)*J11(i, j)   - n(0)*dzdx1(i, j)},
                                                        {n(1)*J12(i, j)   - n(0)*J22(i, j),   n(0)*J21(i, j)   - n(1)*J11(i, j)}};
                        J_red/=sqrt_a(i, j);
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        arma::vec::fixed<2> dmudxi;
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                double  psi2 = w2(j)*t2;
                                double dpsi2 = w2(j)*dt2 + dw2(j)*t2;
                                mu(i, j) += mu_hat(p+q*nx)*t1*psi2;
                                dmudxi += arma::vec::fixed<2>({mu_hat(p+q*nx)*dt1*psi2, mu_hat(p+q*nx)*t1*dpsi2});

                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                        arma::vec::fixed<3> q_mu = J_red*dmudxi;
                        dcp(i, j) = 2*dot(Q, q_mu);
                    }
                break;
            }
            default:
                std::println("Only linear and nonlinear analysis are implemented for Wing!");
                exit(EXIT_FAILURE);
        }
        lift   = 0;
        moment = 0;
        arma::vec  x1_gl = phi1->xg;
        arma::vec  x2_gl = phi2->xg;
        arma::vec  w1_gl = phi1->wg;
        arma::vec  w2_gl = phi2->wg;
        arma::vec dx2_gl = phi2->dxg;
        arma::vec dw2_gl = phi2->dwg;
        arma::vec dw2    = phi2->weightFunctionDerivativeFactor(dx2_gl);
        auto [  x_gl,   y_gl] = Lagrange::TransfiniteQuadMap( x1_gl,  x2_gl, chi);
        auto [d2x_gl, d2y_gl] = Lagrange::TransfiniteQuadMap( x1_gl, dx2_gl, chi);
        auto [ dxdx1_gl,  dxdx2_gl,  dydx1_gl,  dydx2_gl] = Lagrange::TransfiniteQuadMetrics( x1_gl,  x2_gl, chi);
        auto [d2xdx1_gl, d2xdx2_gl, d2ydx1_gl, d2ydx2_gl] = Lagrange::TransfiniteQuadMetrics( x1_gl, dx2_gl, chi);

        switch (analysis)
        {
            case Analysis::linear:
            {
                #pragma omp parallel for reduction(+:lift) reduction(+:moment)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        double detJ  =  dxdx1_gl(i, j)* dydx2_gl(i, j) -  dxdx2_gl(i, j)* dydx1_gl(i, j);
                        double detJ2 = d2xdx1_gl(i, j)*d2ydx2_gl(i, j) - d2xdx2_gl(i, j)*d2ydx1_gl(i, j);
                        double J11_inv  =  dydx2_gl(i, j)/detJ;
                        double J21_inv  = -dydx1_gl(i, j)/detJ;
                        double J21_inv2 =-d2ydx1_gl(i, j)/detJ2;
                        double DCP = 0, DCP2 = 0;
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            double  t2_2   = phi2->constant();
                            double  t2_2p1 = phi2->linear(dx2_gl(j));
                            double dt2_2   = phi2->constantDerivative();
                            double dt2_2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                DCP  += 2*mu_hat(p+q*nx) * (J11_inv*dt1_1*t2_1 + J21_inv*t1_1*dt2_1);
                                DCP2 += 2*mu_hat(p+q*nx) * J21_inv2*dw2(j)*t1_1*t2_2;
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                                std::swap(t2_2, t2_2p1);
                                phi2->next(q, dx2_gl(j), t2_2, t2_2p1);
                                std::swap(dt2_2, dt2_2p1);
                                phi2->nextDerivative(q, dx2_gl(j), t2_2, dt2_2, dt2_2p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                        }
                        double dA  =  w1_gl(i) *  w2_gl(j) * detJ;
                        double dA2 =  w1_gl(i) * dw2_gl(j) * detJ2;
                        lift   += dA * DCP + dA2 * DCP2;
                        moment -= dA * DCP * x_gl(i, j) + dA2 * DCP2 * d2x_gl(i, j);
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec F(3, arma::fill::zeros), M(3, arma::fill::zeros);
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat  Tx_gl = Lagrange::interpolationMatrix(x1,  x1_gl);
                arma::mat  Ty_gl = Lagrange::interpolationMatrix(x2,  x2_gl);
                arma::mat T2y_gl = Lagrange::interpolationMatrix(x2, dx2_gl);
                arma::mat   z_gl = Lagrange::interpolation2D( Tx_gl,  Ty_gl, z,  x1_gl,  x2_gl);
                arma::mat d2z_gl = Lagrange::interpolation2D( Tx_gl, T2y_gl, z,  x1_gl, dx2_gl);
                arma::mat  D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat  D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat Dd2_gl = Lagrange::derivativeMatrix(dx2_gl);
                arma::mat dzdx1_gl = D1_gl*z_gl;
                arma::mat dzdx2_gl = z_gl*D2_gl.t();
                arma::mat d2zdx1_gl = D1_gl*d2z_gl;
                arma::mat d2zdx2_gl = d2z_gl*Dd2_gl.t();
                arma::field<arma::mat> J_gl  = {{ dxdx1_gl,  dxdx2_gl}, { dydx1_gl,  dydx2_gl}};
                arma::field<arma::mat> J2_gl = {{d2xdx1_gl, d2xdx2_gl}, {d2ydx1_gl, d2ydx2_gl}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube e2_c_gl = MetricCo(J2_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::cube ec2_gl  = MetricContra(e2_c_gl);
                arma::mat e_gl  =  e_c_gl.slice(0)%e_c_gl.slice(2)  - pow(e_c_gl.slice(1), 2);
                arma::mat e2_gl = e2_c_gl.slice(0)%e2_c_gl.slice(2) - pow(e2_c_gl.slice(1), 2);
                arma::mat sqrt_a  = sqrt( e_gl%(1 +  ec_gl.slice(0)%pow( dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl   + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
                arma::mat sqrt_a2 = sqrt(e2_gl%(1 + ec2_gl.slice(0)%pow(d2zdx1_gl, 2) + 2*ec2_gl.slice(1)%d2zdx1_gl%d2zdx2_gl + ec2_gl.slice(2)%pow(d2zdx2_gl, 2)));
                #pragma omp parallel for reduction(+:F) reduction(+:M)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        arma::vec::fixed<3> n_gl  = arma::vec::fixed<3>({dydx1_gl(i, j)*dzdx2_gl(i, j)-dzdx1_gl(i, j)*dydx2_gl(i, j),
                                                                         dzdx1_gl(i, j)*dxdx2_gl(i, j)-dxdx1_gl(i, j)*dzdx2_gl(i, j),
                                                                         dxdx1_gl(i, j)*dydx2_gl(i, j)-dydx1_gl(i, j)*dxdx2_gl(i, j)})/sqrt_a(i, j);
                        arma::vec::fixed<3> n2_gl = arma::vec::fixed<3>({d2ydx1_gl(i, j)*d2zdx2_gl(i, j)-d2zdx1_gl(i, j)*d2ydx2_gl(i, j),
                                                                         d2zdx1_gl(i, j)*d2xdx2_gl(i, j)-d2xdx1_gl(i, j)*d2zdx2_gl(i, j),
                                                                         d2xdx1_gl(i, j)*d2ydx2_gl(i, j)-d2ydx1_gl(i, j)*d2xdx2_gl(i, j)})/sqrt_a2(i, j);
                        arma::mat::fixed<3, 2> J_red  = {{dydx2_gl(i, j)*n_gl(2)-n_gl(1)*dzdx2_gl(i, j), n_gl(1)*dzdx1_gl(i, j)-dydx1_gl(i, j)*n_gl(2)},
                                                         {dzdx2_gl(i, j)*n_gl(0)-n_gl(2)*dxdx2_gl(i, j), n_gl(2)*dxdx1_gl(i, j)-dzdx1_gl(i, j)*n_gl(0)},
                                                         {dxdx2_gl(i, j)*n_gl(1)-n_gl(0)*dydx2_gl(i, j), n_gl(0)*dydx1_gl(i, j)-dxdx1_gl(i, j)*n_gl(1)}};
                        arma::mat::fixed<3, 2> J2_red = {{d2ydx2_gl(i, j)*n2_gl(2)-n2_gl(1)*d2zdx2_gl(i, j), n2_gl(1)*d2zdx1_gl(i, j)-d2ydx1_gl(i, j)*n2_gl(2)},
                                                         {d2zdx2_gl(i, j)*n2_gl(0)-n2_gl(2)*d2xdx2_gl(i, j), n2_gl(2)*d2xdx1_gl(i, j)-d2zdx1_gl(i, j)*n2_gl(0)},
                                                         {d2xdx2_gl(i, j)*n2_gl(1)-n2_gl(0)*d2ydx2_gl(i, j), n2_gl(0)*d2ydx1_gl(i, j)-d2xdx1_gl(i, j)*n2_gl(1)}};

                        arma::vec::fixed<2> dmudxi(arma::fill::zeros), dmu1dxi(arma::fill::zeros), dmu2dxi(arma::fill::zeros);
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            double  t2_2   = phi2->constant();
                            double  t2_2p1 = phi2->linear(dx2_gl(j));
                            double dt2_2   = phi2->constantDerivative();
                            double dt2_2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                dmudxi  += mu_hat(p+q*nx)*arma::vec::fixed<2>({dt1_1*t2_1, t1_1*dt2_1});
                                dmu2dxi += mu_hat(p+q*nx)*arma::vec::fixed<2>({0, t1_1*dw2(j)*t2_2});
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                                std::swap(t2_2, t2_2p1);
                                phi2->next(q, dx2_gl(j), t2_2, t2_2p1);
                                std::swap(dt2_2, dt2_2p1);
                                phi2->nextDerivative(q, dx2_gl(j), t2_2, dt2_2, dt2_2p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                        }
                        arma::vec::fixed<3> q_mu  =  J_red*dmudxi;
                        arma::vec::fixed<3> q2_mu = J2_red*dmu2dxi;
                        double DCP  = 2*dot(Q,  q_mu);
                        double DCP2 = 2*dot(Q, q2_mu);
                        arma::vec r  = {  x_gl(i, j),   y_gl(i, j),   z_gl(i, j)};
                        arma::vec r2 = {d2x_gl(i, j), d2y_gl(i, j), d2z_gl(i, j)};
                        F +=  w1_gl(i) *  w2_gl(j) * DCP  *  n_gl
                           +  w1_gl(i) * dw2_gl(j) * DCP2 * n2_gl;
                        M -=  w1_gl(i) * w2_gl(j) * cross( n_gl * DCP,  r)
                            + w1_gl(i) * dw2_gl(j) * cross(n2_gl * DCP2, r2);
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
    else if (phi2->basis == Basis::T)
    {
        mu.zeros();
        dcp.zeros();
        switch (analysis)
        {
            case Analysis::linear:
            {
                arma::mat detJ = J11%J22 - J12%J21;
                arma::mat J11_inv = J22/detJ;
                arma::mat J21_inv =-J21/detJ;
                arma::vec  w1 = phi1->weightFunction(x1);
                arma::vec dw1 = phi1->weightFunctionDerivative(x1);
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  psi1 = w1(i)*t1;
                            double dpsi1 = w1(i)*dt1 + dw1(i)*t1;
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                mu(i, j)  +=   mu_hat(p+q*nx) * psi1*t2;
                                dcp(i, j) += 2*mu_hat(p+q*nx) * (J11_inv(i, j)*dpsi1*t2 + J21_inv(i, j)*psi1*dt2);
                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec  w1 = phi1->weightFunction(x1);
                arma::vec dw1 = phi1->weightFunctionDerivative(x1);
                arma::mat d1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
                arma::mat d2 = Chebyshev::derivativeMatrix(x2, Derivative::first);
                arma::mat dzdx1 = d1 * z;
                arma::mat dzdx2 = z * d2.t();
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat e = e_c.slice(0)%e_c.slice(2) - pow(e_c.slice(1), 2);
                arma::mat sqrt_a = sqrt(e%(1 + ec.slice(0)%pow(dzdx1, 2) + 2*ec.slice(1)%dzdx1%dzdx2 + ec.slice(2)%pow(dzdx2, 2)));
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        arma::vec::fixed<3> n = arma::vec::fixed<3>({J21(i, j)*dzdx2(i, j)-dzdx1(i, j)*J22(i, j),
                                                                     dzdx1(i, j)*J12(i, j)-J11(i, j)*dzdx2(i, j),
                                                                     J11(i, j)*J22(i, j)-J21(i, j)*J12(i, j)})/sqrt_a(i, j);
                        arma::mat::fixed<3, 2> J_red = {{n(2)*J22(i, j)   - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - n(2)*J21(i, j)},
                                                        {n(0)*dzdx2(i, j) - n(2)*J12(i, j),   n(2)*J11(i, j)   - n(0)*dzdx1(i, j)},
                                                        {n(1)*J12(i, j)   - n(0)*J22(i, j),   n(0)*J21(i, j)   - n(1)*J11(i, j)}};
                        J_red/=sqrt_a(i, j);
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        arma::vec::fixed<2> dmudxi;
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  psi1 = w1(i)*t1;
                            double dpsi1 = w1(i)*dt1 + dw1(i)*t1;
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                mu(i, j) += mu_hat(p+q*nx)*psi1*t2;
                                dmudxi += arma::vec::fixed<2>({mu_hat(p+q*nx)*dpsi1*t2, mu_hat(p+q*nx)*psi1*dt2});

                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                        arma::vec::fixed<3> q_mu = J_red*dmudxi;
                        dcp(i, j) = 2*dot(Q, q_mu);
                    }
                break;
            }
            default:
                std::println("Only linear and nonlinear analysis are implemented for Wing!");
                exit(EXIT_FAILURE);
        }
        lift   = 0;
        moment = 0;
        arma::vec  x1_gl = phi1->xg;
        arma::vec  x2_gl = phi2->xg;
        arma::vec  w1_gl = phi1->wg;
        arma::vec  w2_gl = phi2->wg;
        arma::vec dx1_gl = phi1->dxg;
        arma::vec dw1_gl = phi1->dwg;
        arma::vec dw1    = phi1->weightFunctionDerivativeFactor(dx1_gl);
        auto [  x_gl,   y_gl] = Lagrange::TransfiniteQuadMap( x1_gl,  x2_gl, chi);
        auto [d1x_gl, d1y_gl] = Lagrange::TransfiniteQuadMap(dx1_gl,  x2_gl, chi);
        auto [ dxdx1_gl,  dxdx2_gl,  dydx1_gl,  dydx2_gl] = Lagrange::TransfiniteQuadMetrics( x1_gl,  x2_gl, chi);
        auto [d1xdx1_gl, d1xdx2_gl, d1ydx1_gl, d1ydx2_gl] = Lagrange::TransfiniteQuadMetrics(dx1_gl,  x2_gl, chi);

        switch (analysis)
        {
            case Analysis::linear:
            {
                #pragma omp parallel for reduction(+:lift) reduction(+:moment)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        double detJ  =  dxdx1_gl(i, j)* dydx2_gl(i, j) -  dxdx2_gl(i, j)* dydx1_gl(i, j);
                        double detJ1 = d1xdx1_gl(i, j)*d1ydx2_gl(i, j) - d1xdx2_gl(i, j)*d1ydx1_gl(i, j);
                        double J11_inv  =  dydx2_gl(i, j)/detJ;
                        double J21_inv  = -dydx1_gl(i, j)/detJ;
                        double J11_inv1 = d1ydx2_gl(i, j)/detJ1;
                        double DCP = 0, DCP1 = 0;
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        double  t1_2   = phi1->constant();
                        double  t1_2p1 = phi1->linear(dx1_gl(i));
                        double dt1_2   = phi1->constantDerivative();
                        double dt1_2p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                DCP  += 2*mu_hat(p+q*nx) * (J11_inv*dt1_1*t2_1 + J21_inv*t1_1*dt2_1);
                                DCP1 += 2*mu_hat(p+q*nx) * J11_inv1*dw1(i)*t1_2*t2_1;
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                            std::swap(t1_2, t1_2p1);
                            phi1->next(p, dx1_gl(i), t1_2, t1_2p1);
                            std::swap(dt1_2, dt1_2p1);
                            phi1->nextDerivative(p, dx1_gl(i), t1_2, dt1_2, dt1_2p1);
                        }
                        double dA  =  w1_gl(i) *  w2_gl(j) * detJ;
                        double dA1 = dw1_gl(i) *  w2_gl(j) * detJ1;
                        lift   += dA * DCP + dA1 * DCP1;
                        moment -= dA * DCP * x_gl(i, j) + dA1 * DCP1 * d1x_gl(i, j);
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec F(3, arma::fill::zeros), M(3, arma::fill::zeros);
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat  Tx_gl = Lagrange::interpolationMatrix(x1,  x1_gl);
                arma::mat  Ty_gl = Lagrange::interpolationMatrix(x2,  x2_gl);
                arma::mat T1x_gl = Lagrange::interpolationMatrix(x1, dx1_gl);
                arma::mat   z_gl = Lagrange::interpolation2D( Tx_gl,  Ty_gl, z,  x1_gl,  x2_gl);
                arma::mat d1z_gl = Lagrange::interpolation2D(T1x_gl,  Ty_gl, z, dx1_gl,  x2_gl);
                arma::mat  D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat  D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat Dd1_gl = Lagrange::derivativeMatrix(dx1_gl);
                arma::mat dzdx1_gl = D1_gl*z_gl;
                arma::mat dzdx2_gl = z_gl*D2_gl.t();
                arma::mat d1zdx1_gl = Dd1_gl*d1z_gl;
                arma::mat d1zdx2_gl = d1z_gl*D2_gl.t();
                arma::field<arma::mat> J_gl  = {{ dxdx1_gl,  dxdx2_gl}, { dydx1_gl,  dydx2_gl}};
                arma::field<arma::mat> J1_gl = {{d1xdx1_gl, d1xdx2_gl}, {d1ydx1_gl, d1ydx2_gl}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube e1_c_gl = MetricCo(J1_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::cube ec1_gl  = MetricContra(e1_c_gl);
                arma::mat e_gl  =  e_c_gl.slice(0)%e_c_gl.slice(2)  - pow(e_c_gl.slice(1), 2);
                arma::mat e1_gl = e1_c_gl.slice(0)%e1_c_gl.slice(2) - pow(e1_c_gl.slice(1), 2);
                arma::mat sqrt_a  = sqrt( e_gl%(1 +  ec_gl.slice(0)%pow( dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl   + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
                arma::mat sqrt_a1 = sqrt(e1_gl%(1 + ec1_gl.slice(0)%pow(d1zdx1_gl, 2) + 2*ec1_gl.slice(1)%d1zdx1_gl%d1zdx2_gl + ec1_gl.slice(2)%pow(d1zdx2_gl, 2)));
                #pragma omp parallel for reduction(+:F) reduction(+:M)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        arma::vec::fixed<3> n_gl  = arma::vec::fixed<3>({dydx1_gl(i, j)*dzdx2_gl(i, j)-dzdx1_gl(i, j)*dydx2_gl(i, j),
                                                                         dzdx1_gl(i, j)*dxdx2_gl(i, j)-dxdx1_gl(i, j)*dzdx2_gl(i, j),
                                                                         dxdx1_gl(i, j)*dydx2_gl(i, j)-dydx1_gl(i, j)*dxdx2_gl(i, j)})/sqrt_a(i, j);
                        arma::vec::fixed<3> n1_gl = arma::vec::fixed<3>({d1ydx1_gl(i, j)*d1zdx2_gl(i, j)-d1zdx1_gl(i, j)*d1ydx2_gl(i, j),
                                                                         d1zdx1_gl(i, j)*d1xdx2_gl(i, j)-d1xdx1_gl(i, j)*d1zdx2_gl(i, j),
                                                                         d1xdx1_gl(i, j)*d1ydx2_gl(i, j)-d1ydx1_gl(i, j)*d1xdx2_gl(i, j)})/sqrt_a1(i, j);
                        arma::mat::fixed<3, 2> J_red  = {{dydx2_gl(i, j)*n_gl(2)-n_gl(1)*dzdx2_gl(i, j), n_gl(1)*dzdx1_gl(i, j)-dydx1_gl(i, j)*n_gl(2)},
                                                         {dzdx2_gl(i, j)*n_gl(0)-n_gl(2)*dxdx2_gl(i, j), n_gl(2)*dxdx1_gl(i, j)-dzdx1_gl(i, j)*n_gl(0)},
                                                         {dxdx2_gl(i, j)*n_gl(1)-n_gl(0)*dydx2_gl(i, j), n_gl(0)*dydx1_gl(i, j)-dxdx1_gl(i, j)*n_gl(1)}};
                        arma::mat::fixed<3, 2> J1_red = {{d1ydx2_gl(i, j)*n1_gl(2)-n1_gl(1)*d1zdx2_gl(i, j), n1_gl(1)*d1zdx1_gl(i, j)-d1ydx1_gl(i, j)*n1_gl(2)},
                                                         {d1zdx2_gl(i, j)*n1_gl(0)-n1_gl(2)*d1xdx2_gl(i, j), n1_gl(2)*d1xdx1_gl(i, j)-d1zdx1_gl(i, j)*n1_gl(0)},
                                                         {d1xdx2_gl(i, j)*n1_gl(1)-n1_gl(0)*d1ydx2_gl(i, j), n1_gl(0)*d1ydx1_gl(i, j)-d1xdx1_gl(i, j)*n1_gl(1)}};

                        arma::vec::fixed<2> dmudxi(arma::fill::zeros), dmu1dxi(arma::fill::zeros), dmu2dxi(arma::fill::zeros);
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        double  t1_2   = phi1->constant();
                        double  t1_2p1 = phi1->linear(dx1_gl(i));
                        double dt1_2   = phi1->constantDerivative();
                        double dt1_2p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                dmudxi  += mu_hat(p+q*nx)*arma::vec::fixed<2>({dt1_1*t2_1, t1_1*dt2_1});
                                dmu1dxi += mu_hat(p+q*nx)*arma::vec::fixed<2>({dw1(i)*t1_2*t2_1, 0});
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                            std::swap(t1_2, t1_2p1);
                            phi1->next(p, dx1_gl(i), t1_2, t1_2p1);
                            std::swap(dt1_2, dt1_2p1);
                            phi1->nextDerivative(p, dx1_gl(i), t1_2, dt1_2, dt1_2p1);
                        }
                        arma::vec::fixed<3> q_mu  =  J_red*dmudxi;
                        arma::vec::fixed<3> q1_mu = J1_red*dmu1dxi;
                        double DCP  = 2*dot(Q,  q_mu);
                        double DCP1 = 2*dot(Q, q1_mu);
                        arma::vec r  = {  x_gl(i, j),   y_gl(i, j),   z_gl(i, j)};
                        arma::vec r1 = {d1x_gl(i, j), d1y_gl(i, j), d1z_gl(i, j)};
                        F +=  w1_gl(i) * w2_gl(j) * DCP  *  n_gl
                           + dw1_gl(i) * w2_gl(j) * DCP1 * n1_gl;
                        M -=  w1_gl(i) * w2_gl(j) * cross( n_gl * DCP,  r)
                           + dw1_gl(i) * w2_gl(j) * cross(n1_gl * DCP1, r1);
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
    else
    {
        mu.zeros();
        dcp.zeros();
        switch (analysis)
        {
            case Analysis::linear:
            {
                arma::mat detJ = J11%J22 - J12%J21;
                arma::mat J11_inv = J22/detJ;
                arma::mat J21_inv =-J21/detJ;
                arma::vec  w1 = phi1->weightFunction(x1);
                arma::vec  w2 = phi2->weightFunction(x2);
                arma::vec dw1 = phi1->weightFunctionDerivative(x1);
                arma::vec dw2 = phi2->weightFunctionDerivative(x2);
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  psi1 = w1(i)*t1;
                            double dpsi1 = w1(i)*dt1 + dw1(i)*t1;
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                double  psi2 = w2(j)*t2;
                                double dpsi2 = w2(j)*dt2 + dw2(j)*t2;
                                mu(i, j)  +=   mu_hat(p+q*nx) * psi1*psi2;
                                dcp(i, j) += 2*mu_hat(p+q*nx) * (J11_inv(i, j)*dpsi1*psi2 + J21_inv(i, j)*psi1*dpsi2);
                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec  w1 = phi1->weightFunction(x1);
                arma::vec  w2 = phi2->weightFunction(x2);
                arma::vec dw1 = phi1->weightFunctionDerivative(x1);
                arma::vec dw2 = phi2->weightFunctionDerivative(x2);
                arma::mat d1 = Chebyshev::derivativeMatrix(x1, Derivative::first);
                arma::mat d2 = Chebyshev::derivativeMatrix(x2, Derivative::first);
                arma::mat dzdx1 = d1 * z;
                arma::mat dzdx2 = z * d2.t();
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat e = e_c.slice(0)%e_c.slice(2) - pow(e_c.slice(1), 2);
                arma::mat sqrt_a = sqrt(e%(1 + ec.slice(0)%pow(dzdx1, 2) + 2*ec.slice(1)%dzdx1%dzdx2 + ec.slice(2)%pow(dzdx2, 2)));
                #pragma omp parallel for
                for (size_t i = 0; i < nx; i++) // Loop over nodes in 1-direction
                    for (size_t j = 0; j < ny; j++) // Loop over nodes in 2-direction
                    {
                        arma::vec::fixed<3> n = arma::vec::fixed<3>({J21(i, j)*dzdx2(i, j) - J22(i, j)*dzdx1(i, j),
                                                                     J12(i, j)*dzdx1(i, j) - J11(i, j)*dzdx2(i, j),
                                                                     J11(i, j)*J22(i, j)   - J21(i, j)*J12(i, j)})/sqrt_a(i, j);
                        arma::mat::fixed<3, 2> J_red = {{n(2)*J22(i, j)   - n(1)*dzdx2(i, j), n(1)*dzdx1(i, j) - n(2)*J21(i, j)},
                                                        {n(0)*dzdx2(i, j) - n(2)*J12(i, j),   n(2)*J11(i, j)   - n(0)*dzdx1(i, j)},
                                                        {n(1)*J12(i, j)   - n(0)*J22(i, j),   n(0)*J21(i, j)   - n(1)*J11(i, j)}};
                        J_red/=sqrt_a(i, j);
                        double  t1   = phi1->constant();
                        double  t1p1 = phi1->linear(x1(i));
                        double dt1   = phi1->constantDerivative();
                        double dt1p1 = phi1->linearDerivative();
                        arma::vec::fixed<2> dmudxi;
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  psi1 = w1(i)*t1;
                            double dpsi1 = w1(i)*dt1 + dw1(i)*t1;
                            double  t2   = phi2->constant();
                            double  t2p1 = phi2->linear(x2(j));
                            double dt2   = phi2->constantDerivative();
                            double dt2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                double  psi2 = w2(j)*t2;
                                double dpsi2 = w2(j)*dt2 + dw2(j)*t2;
                                mu(i, j) += mu_hat(p+q*nx)*psi1*psi2;
                                dmudxi += arma::vec::fixed<2>({mu_hat(p+q*nx)*dpsi1*psi2, mu_hat(p+q*nx)*psi1*dpsi2});

                                std::swap(t2, t2p1);
                                phi2->next(q, x2(j), t2, t2p1);
                                std::swap(dt2, dt2p1);
                                phi2->nextDerivative(q, x2(j), t2, dt2, dt2p1);
                            }
                            std::swap(t1, t1p1);
                            phi1->next(p, x1(i), t1, t1p1);
                            std::swap(dt1, dt1p1);
                            phi1->nextDerivative(p, x1(i), t1, dt1, dt1p1);
                        }
                        arma::vec::fixed<3> q_mu = J_red*dmudxi;
                        dcp(i, j) = 2*dot(Q, q_mu);
                    }
                break;
            }
            default:
                std::println("Only linear and nonlinear analysis are implemented for Wing!");
                exit(EXIT_FAILURE);
        }
        lift   = 0;
        moment = 0;
        arma::vec  x1_gl = phi1->xg;
        arma::vec  x2_gl = phi2->xg;
        arma::vec  w1_gl = phi1->wg;
        arma::vec  w2_gl = phi2->wg;
        arma::vec dx1_gl = phi1->dxg;
        arma::vec dx2_gl = phi2->dxg;
        arma::vec dw1_gl = phi1->dwg;
        arma::vec dw2_gl = phi2->dwg;
        arma::vec dw1    = phi1->weightFunctionDerivativeFactor(dx1_gl);
        arma::vec dw2    = phi2->weightFunctionDerivativeFactor(dx2_gl);
        auto [  x_gl,   y_gl] = Lagrange::TransfiniteQuadMap( x1_gl,  x2_gl, chi);
        auto [d1x_gl, d1y_gl] = Lagrange::TransfiniteQuadMap(dx1_gl,  x2_gl, chi);
        auto [d2x_gl, d2y_gl] = Lagrange::TransfiniteQuadMap( x1_gl, dx2_gl, chi);
        auto [ dxdx1_gl,  dxdx2_gl,  dydx1_gl,  dydx2_gl] = Lagrange::TransfiniteQuadMetrics( x1_gl,  x2_gl, chi);
        auto [d1xdx1_gl, d1xdx2_gl, d1ydx1_gl, d1ydx2_gl] = Lagrange::TransfiniteQuadMetrics(dx1_gl,  x2_gl, chi);
        auto [d2xdx1_gl, d2xdx2_gl, d2ydx1_gl, d2ydx2_gl] = Lagrange::TransfiniteQuadMetrics( x1_gl, dx2_gl, chi);

        switch (analysis)
        {
            case Analysis::linear:
            {
                #pragma omp parallel for reduction(+:lift) reduction(+:moment)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        double detJ  =  dxdx1_gl(i, j)* dydx2_gl(i, j) -  dxdx2_gl(i, j)* dydx1_gl(i, j);
                        double detJ1 = d1xdx1_gl(i, j)*d1ydx2_gl(i, j) - d1xdx2_gl(i, j)*d1ydx1_gl(i, j);
                        double detJ2 = d2xdx1_gl(i, j)*d2ydx2_gl(i, j) - d2xdx2_gl(i, j)*d2ydx1_gl(i, j);
                        double J11_inv  =  dydx2_gl(i, j)/detJ;
                        double J21_inv  = -dydx1_gl(i, j)/detJ;
                        double J11_inv1 = d1ydx2_gl(i, j)/detJ1;
                        double J21_inv2 =-d2ydx1_gl(i, j)/detJ2;
                        double DCP = 0, DCP1 = 0, DCP2 = 0;
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        double  t1_2   = phi1->constant();
                        double  t1_2p1 = phi1->linear(dx1_gl(i));
                        double dt1_2   = phi1->constantDerivative();
                        double dt1_2p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            double  t2_2   = phi2->constant();
                            double  t2_2p1 = phi2->linear(dx2_gl(j));
                            double dt2_2   = phi2->constantDerivative();
                            double dt2_2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                DCP  += 2*mu_hat(p+q*nx) * (J11_inv*dt1_1*t2_1 + J21_inv*t1_1*dt2_1);
                                DCP1 += 2*mu_hat(p+q*nx) * J11_inv1*dw1(i)*t1_2*t2_1;
                                DCP2 += 2*mu_hat(p+q*nx) * J21_inv2*dw2(j)*t1_1*t2_2;
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                                std::swap(t2_2, t2_2p1);
                                phi2->next(q, dx2_gl(j), t2_2, t2_2p1);
                                std::swap(dt2_2, dt2_2p1);
                                phi2->nextDerivative(q, dx2_gl(j), t2_2, dt2_2, dt2_2p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                            std::swap(t1_2, t1_2p1);
                            phi1->next(p, dx1_gl(i), t1_2, t1_2p1);
                            std::swap(dt1_2, dt1_2p1);
                            phi1->nextDerivative(p, dx1_gl(i), t1_2, dt1_2, dt1_2p1);
                        }
                        double dA  =  w1_gl(i) *  w2_gl(j) * detJ;
                        double dA1 = dw1_gl(i) *  w2_gl(j) * detJ1;
                        double dA2 =  w1_gl(i) * dw2_gl(j) * detJ2;
                        lift   += dA * DCP + dA1 * DCP1 + dA2 * DCP2;
                        moment -= dA * DCP * x_gl(i, j) + dA1 * DCP1 * d1x_gl(i, j) + dA2 * DCP2 * d2x_gl(i, j);
                    }
                break;
            }
            case Analysis::nonlinear:
            {
                arma::vec F(3, arma::fill::zeros), M(3, arma::fill::zeros);
                arma::vec Q = {cos(alpha), 0, sin(alpha)};
                arma::mat  Tx_gl = Lagrange::interpolationMatrix(x1,  x1_gl);
                arma::mat  Ty_gl = Lagrange::interpolationMatrix(x2,  x2_gl);
                arma::mat T1x_gl = Lagrange::interpolationMatrix(x1, dx1_gl);
                arma::mat T2y_gl = Lagrange::interpolationMatrix(x2, dx2_gl);
                arma::mat   z_gl = Lagrange::interpolation2D( Tx_gl,  Ty_gl, z,  x1_gl,  x2_gl);
                arma::mat d1z_gl = Lagrange::interpolation2D(T1x_gl,  Ty_gl, z, dx1_gl,  x2_gl);
                arma::mat d2z_gl = Lagrange::interpolation2D( Tx_gl, T2y_gl, z,  x1_gl, dx2_gl);
                arma::mat  D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat  D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat Dd1_gl = Lagrange::derivativeMatrix(dx1_gl);
                arma::mat Dd2_gl = Lagrange::derivativeMatrix(dx2_gl);
                arma::mat dzdx1_gl = D1_gl*z_gl;
                arma::mat dzdx2_gl = z_gl*D2_gl.t();
                arma::mat d1zdx1_gl = Dd1_gl*d1z_gl;
                arma::mat d2zdx1_gl = D1_gl*d2z_gl;
                arma::mat d1zdx2_gl = d1z_gl*D2_gl.t();
                arma::mat d2zdx2_gl = d2z_gl*Dd2_gl.t();
                arma::field<arma::mat> J_gl  = {{ dxdx1_gl,  dxdx2_gl}, { dydx1_gl,  dydx2_gl}};
                arma::field<arma::mat> J1_gl = {{d1xdx1_gl, d1xdx2_gl}, {d1ydx1_gl, d1ydx2_gl}};
                arma::field<arma::mat> J2_gl = {{d2xdx1_gl, d2xdx2_gl}, {d2ydx1_gl, d2ydx2_gl}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube e1_c_gl = MetricCo(J1_gl);
                arma::cube e2_c_gl = MetricCo(J2_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::cube ec1_gl  = MetricContra(e1_c_gl);
                arma::cube ec2_gl  = MetricContra(e2_c_gl);
                arma::mat e_gl  =  e_c_gl.slice(0)%e_c_gl.slice(2)  - pow(e_c_gl.slice(1), 2);
                arma::mat e1_gl = e1_c_gl.slice(0)%e1_c_gl.slice(2) - pow(e1_c_gl.slice(1), 2);
                arma::mat e2_gl = e2_c_gl.slice(0)%e2_c_gl.slice(2) - pow(e2_c_gl.slice(1), 2);
                arma::mat sqrt_a  = sqrt( e_gl%(1 +  ec_gl.slice(0)%pow( dzdx1_gl, 2) + 2*ec_gl.slice(1)%dzdx1_gl%dzdx2_gl   + ec_gl.slice(2)%pow(dzdx2_gl, 2)));
                arma::mat sqrt_a1 = sqrt(e1_gl%(1 + ec1_gl.slice(0)%pow(d1zdx1_gl, 2) + 2*ec1_gl.slice(1)%d1zdx1_gl%d1zdx2_gl + ec1_gl.slice(2)%pow(d1zdx2_gl, 2)));
                arma::mat sqrt_a2 = sqrt(e2_gl%(1 + ec2_gl.slice(0)%pow(d2zdx1_gl, 2) + 2*ec2_gl.slice(1)%d2zdx1_gl%d2zdx2_gl + ec2_gl.slice(2)%pow(d2zdx2_gl, 2)));
                #pragma omp parallel for reduction(+:F) reduction(+:M)
                for (size_t i = 0; i < nx; i++)
                    for (size_t j = 0; j < ny; j++)
                    {
                        arma::vec::fixed<3> n_gl  = arma::vec::fixed<3>({dydx1_gl(i, j)*dzdx2_gl(i, j)-dzdx1_gl(i, j)*dydx2_gl(i, j),
                                                                         dzdx1_gl(i, j)*dxdx2_gl(i, j)-dxdx1_gl(i, j)*dzdx2_gl(i, j),
                                                                         dxdx1_gl(i, j)*dydx2_gl(i, j)-dydx1_gl(i, j)*dxdx2_gl(i, j)})/sqrt_a(i, j);
                        arma::vec::fixed<3> n1_gl = arma::vec::fixed<3>({d1ydx1_gl(i, j)*d1zdx2_gl(i, j)-d1zdx1_gl(i, j)*d1ydx2_gl(i, j),
                                                                         d1zdx1_gl(i, j)*d1xdx2_gl(i, j)-d1xdx1_gl(i, j)*d1zdx2_gl(i, j),
                                                                         d1xdx1_gl(i, j)*d1ydx2_gl(i, j)-d1ydx1_gl(i, j)*d1xdx2_gl(i, j)})/sqrt_a1(i, j);
                        arma::vec::fixed<3> n2_gl = arma::vec::fixed<3>({d2ydx1_gl(i, j)*d2zdx2_gl(i, j)-d2zdx1_gl(i, j)*d2ydx2_gl(i, j),
                                                                         d2zdx1_gl(i, j)*d2xdx2_gl(i, j)-d2xdx1_gl(i, j)*d2zdx2_gl(i, j),
                                                                         d2xdx1_gl(i, j)*d2ydx2_gl(i, j)-d2ydx1_gl(i, j)*d2xdx2_gl(i, j)})/sqrt_a2(i, j);
                        arma::mat::fixed<3, 2> J_red  = {{dydx2_gl(i, j)*n_gl(2)-n_gl(1)*dzdx2_gl(i, j), n_gl(1)*dzdx1_gl(i, j)-dydx1_gl(i, j)*n_gl(2)},
                                                         {dzdx2_gl(i, j)*n_gl(0)-n_gl(2)*dxdx2_gl(i, j), n_gl(2)*dxdx1_gl(i, j)-dzdx1_gl(i, j)*n_gl(0)},
                                                         {dxdx2_gl(i, j)*n_gl(1)-n_gl(0)*dydx2_gl(i, j), n_gl(0)*dydx1_gl(i, j)-dxdx1_gl(i, j)*n_gl(1)}};
                        arma::mat::fixed<3, 2> J1_red = {{d1ydx2_gl(i, j)*n1_gl(2)-n1_gl(1)*d1zdx2_gl(i, j), n1_gl(1)*d1zdx1_gl(i, j)-d1ydx1_gl(i, j)*n1_gl(2)},
                                                         {d1zdx2_gl(i, j)*n1_gl(0)-n1_gl(2)*d1xdx2_gl(i, j), n1_gl(2)*d1xdx1_gl(i, j)-d1zdx1_gl(i, j)*n1_gl(0)},
                                                         {d1xdx2_gl(i, j)*n1_gl(1)-n1_gl(0)*d1ydx2_gl(i, j), n1_gl(0)*d1ydx1_gl(i, j)-d1xdx1_gl(i, j)*n1_gl(1)}};
                        arma::mat::fixed<3, 2> J2_red = {{d2ydx2_gl(i, j)*n2_gl(2)-n2_gl(1)*d2zdx2_gl(i, j), n2_gl(1)*d2zdx1_gl(i, j)-d2ydx1_gl(i, j)*n2_gl(2)},
                                                         {d2zdx2_gl(i, j)*n2_gl(0)-n2_gl(2)*d2xdx2_gl(i, j), n2_gl(2)*d2xdx1_gl(i, j)-d2zdx1_gl(i, j)*n2_gl(0)},
                                                         {d2xdx2_gl(i, j)*n2_gl(1)-n2_gl(0)*d2ydx2_gl(i, j), n2_gl(0)*d2ydx1_gl(i, j)-d2xdx1_gl(i, j)*n2_gl(1)}};

                        arma::vec::fixed<2> dmudxi(arma::fill::zeros), dmu1dxi(arma::fill::zeros), dmu2dxi(arma::fill::zeros);
                        double  t1_1   = phi1->constant();
                        double  t1_1p1 = phi1->linear(x1_gl(i));
                        double dt1_1   = phi1->constantDerivative();
                        double dt1_1p1 = phi1->linearDerivative();
                        double  t1_2   = phi1->constant();
                        double  t1_2p1 = phi1->linear(dx1_gl(i));
                        double dt1_2   = phi1->constantDerivative();
                        double dt1_2p1 = phi1->linearDerivative();
                        for (size_t p = 0; p < nx; p++)
                        {
                            double  t2_1   = phi2->constant();
                            double  t2_1p1 = phi2->linear(x2_gl(j));
                            double dt2_1   = phi2->constantDerivative();
                            double dt2_1p1 = phi2->linearDerivative();
                            double  t2_2   = phi2->constant();
                            double  t2_2p1 = phi2->linear(dx2_gl(j));
                            double dt2_2   = phi2->constantDerivative();
                            double dt2_2p1 = phi2->linearDerivative();
                            for (size_t q = 0; q < ny; q++)
                            {
                                dmudxi  += mu_hat(p+q*nx)*arma::vec::fixed<2>({dt1_1*t2_1, t1_1*dt2_1});
                                dmu1dxi += mu_hat(p+q*nx)*arma::vec::fixed<2>({dw1(i)*t1_2*t2_1, 0});
                                dmu2dxi += mu_hat(p+q*nx)*arma::vec::fixed<2>({0, t1_1*dw2(j)*t2_2});
                                std::swap(t2_1, t2_1p1);
                                phi2->next(q, x2_gl(j), t2_1, t2_1p1);
                                std::swap(dt2_1, dt2_1p1);
                                phi2->nextDerivative(q, x2_gl(j), t2_1, dt2_1, dt2_1p1);
                                std::swap(t2_2, t2_2p1);
                                phi2->next(q, dx2_gl(j), t2_2, t2_2p1);
                                std::swap(dt2_2, dt2_2p1);
                                phi2->nextDerivative(q, dx2_gl(j), t2_2, dt2_2, dt2_2p1);
                            }
                            std::swap(t1_1, t1_1p1);
                            phi1->next(p, x1_gl(i), t1_1, t1_1p1);
                            std::swap(dt1_1, dt1_1p1);
                            phi1->nextDerivative(p, x1_gl(i), t1_1, dt1_1, dt1_1p1);
                            std::swap(t1_2, t1_2p1);
                            phi1->next(p, dx1_gl(i), t1_2, t1_2p1);
                            std::swap(dt1_2, dt1_2p1);
                            phi1->nextDerivative(p, dx1_gl(i), t1_2, dt1_2, dt1_2p1);
                        }
                        arma::vec::fixed<3> q_mu  =  J_red*dmudxi;
                        arma::vec::fixed<3> q1_mu = J1_red*dmu1dxi;
                        arma::vec::fixed<3> q2_mu = J2_red*dmu2dxi;
                        double DCP  = 2*dot(Q,  q_mu);
                        double DCP1 = 2*dot(Q, q1_mu);
                        double DCP2 = 2*dot(Q, q2_mu);
                        arma::vec r  = {  x_gl(i, j),   y_gl(i, j),   z_gl(i, j)};
                        arma::vec r1 = {d1x_gl(i, j), d1y_gl(i, j), d1z_gl(i, j)};
                        arma::vec r2 = {d2x_gl(i, j), d2y_gl(i, j), d2z_gl(i, j)};
                        F +=  w1_gl(i) *  w2_gl(j) * DCP  *  n_gl
                           + dw1_gl(i) *  w2_gl(j) * DCP1 * n1_gl
                           +  w1_gl(i) * dw2_gl(j) * DCP2 * n2_gl;
                        M -=  w1_gl(i) *  w2_gl(j) * cross( n_gl * DCP,  r)
                           + dw1_gl(i) *  w2_gl(j) * cross(n1_gl * DCP1, r1)
                           +  w1_gl(i) * dw2_gl(j) * cross(n2_gl * DCP2, r2);
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
}