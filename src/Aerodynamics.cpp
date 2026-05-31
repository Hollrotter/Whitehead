#include "Aerodynamics.hpp"

Aerodynamics Aerodynamics::fromWings(std::vector<Wing*> _wings)
{
    std::vector<Interface> _interfaces;
    for (size_t sD = 0; sD < _wings.size()-1; sD++)
        for (size_t sC = 0; sC < 4; sC++)
            for (size_t tD = sD+1; tD < _wings.size(); tD++)
                for (size_t tC = 0; tC < 4 ; tC++)
                    if (_wings[sD]->chi[sC] == _wings[tD]->chi[tC])
                    {
                        _interfaces.push_back(Interface(sD, tD, sC, tC));
                        _wings[sD]->chi[sC]->curveType = CurveType::Interface;
                    }
    double lambda0 = 2;
    for (Interface& interface:_interfaces)
    {
        switch (interface.sourceCurve)
        {
            case 0: // South
                interface.lambdaSource = lambda0*mean(_wings[interface.sourceDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaSource = lambda0*mean(_wings[interface.sourceDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaSource = lambda0*mean(_wings[interface.sourceDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaSource = lambda0*mean(_wings[interface.sourceDomain]->h_1s1_west);
                break;
        }
        switch (interface.targetCurve)
        {
            case 0: // South
                interface.lambdaTarget = lambda0*mean(_wings[interface.targetDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaTarget = lambda0*mean(_wings[interface.targetDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaTarget = lambda0*mean(_wings[interface.targetDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaTarget = lambda0*mean(_wings[interface.targetDomain]->h_1s1_west);
                break;
        }
    }
    return {_wings, _interfaces};
}

void Aerodynamics::setlambda(double l)
{
    lambda0 = l;
    for (Interface& interface:interfaces)
    {
        switch (interface.sourceCurve)
        {
            case 0: // South
                interface.lambdaSource = lambda0*mean(wings[interface.sourceDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaSource = lambda0*mean(wings[interface.sourceDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaSource = lambda0*mean(wings[interface.sourceDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaSource = lambda0*mean(wings[interface.sourceDomain]->h_1s1_west);
                break;
        }
        switch (interface.targetCurve)
        {
            case 0: // South
                interface.lambdaTarget = lambda0*mean(wings[interface.targetDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaTarget = lambda0*mean(wings[interface.targetDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaTarget = lambda0*mean(wings[interface.targetDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaTarget = lambda0*mean(wings[interface.targetDomain]->h_1s1_west);
                break;
        }
    }
}

void Aerodynamics::checkMesh()
{
    for (size_t i = 0; i < wings.size(); i++)
    {
        printf("Checking mesh of wing number %lu\n", i);
        wings[i]->checkMesh();
    }
}

template <class C> void Aerodynamics::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const C val)
{
    for (auto &wing:wings)
        for (auto &chi:wing->chi)
            if (chi == dir)
                wing->boundary(chi, bc, val);
}

template <class C> void Aerodynamics::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2, const C val)
{
    for (auto &wing:wings)
        for (auto &chi:wing->chi)
            if (chi == dir)
                wing->boundary(chi, bc, _r1, _r2, val);
}

void Aerodynamics::linear()
{
    analysis = Analysis::linear;
    // Influence of the wing surfaces on each other
    bw.set_size(wings.size(), wings.size());
    for (size_t sD = 0; sD < wings.size(); sD++)
    {
        wings[sD]->analysis = Analysis::linear;
        for (size_t tD = 0; tD < wings.size(); tD++)
            if (sD != tD)
            {
                std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                {
                    gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                    x1_gl(ii) =-gl_x[ii].x();
                }
                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                {
                    gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                    x2_gl(jj) =-gl_y[jj].x();
                }
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                arma::mat xC = wings[tD]->xC;
                arma::mat yC = wings[tD]->yC;
                bw(tD, sD).set_size(wings[tD]->nxy, wings[sD]->nxy, wings[sD]->con);
                for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                    for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                    {
                        size_t k = i+j*wings[tD]->nx;
                        for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                            for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                            {
                                arma::vec::fixed<2> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j)};
                                double r3 = pow(norm(r), 3);
                                for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction of source
                                {
                                    double t2  = boost::math::chebyshev_t(q, x2_gl(jj));
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                    for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction of source
                                    {
                                        double t1  = boost::math::chebyshev_t(p, x1_gl(ii));
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                        double dmudx1 = dt1 *  t2;
                                        double dmudx2 =  t1 * dt2;
                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= gl_x[ii].weight * gl_y[jj].weight / r3
                                            *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                            - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                                    }
                                }
                            }
                    }
            }
    }
    if (sym == Symmetry::y)
    {
        for (size_t sD = 0; sD < wings.size(); sD++)
            for (size_t tD = 0; tD < wings.size(); tD++)
                if (sD != tD)
                {
                    std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                    arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                    {
                        gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                        x1_gl(ii) =-gl_x[ii].x();
                    }
                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                    {
                        gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                        x2_gl(jj) =-gl_y[jj].x();
                    }
                    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                    y_gl     =-y_gl;
                    dy_gldx1 =-dy_gldx1;
                    dy_gldx2 =-dy_gldx2;
                    arma::mat xC = wings[tD]->xC;
                    arma::mat yC = wings[tD]->yC;
                    for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                        for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                        {
                            size_t k = i+j*wings[tD]->nx;
                            for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                                {
                                    arma::vec::fixed<2> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j)};
                                    double r3 = pow(norm(r), 3);
                                    for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction of source
                                    {
                                        double t2  = boost::math::chebyshev_t(q, x2_gl(jj));
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                        for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction of source
                                        {
                                            double t1  = boost::math::chebyshev_t(p, x1_gl(ii));
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                            double dmudx1 = dt1 *  t2;
                                            double dmudx2 =  t1 * dt2;
                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= gl_x[ii].weight * gl_y[jj].weight / r3
                                                *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                                - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                                        }
                                    }
                                }
                        }
                }
    }
    solve();
}

void Aerodynamics::nonlinear()
{
    analysis = Analysis::nonlinear;
    // Influence of the wing surfaces on each other
    bw.set_size(wings.size(), wings.size());
    for (size_t sD = 0; sD < wings.size(); sD++)
    {
        wings[sD]->analysis = Analysis::nonlinear;
        for (size_t tD = 0; tD < wings.size(); tD++)
            if (sD != tD)
            {
                std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                {
                    gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                    x1_gl(ii) =-gl_x[ii].x();
                }
                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                {
                    gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                    x2_gl(jj) =-gl_y[jj].x();
                }
                arma::mat Tx = Lagrange::interpolationMatrix(wings[sD]->x1, x1_gl);
                arma::mat Ty = Lagrange::interpolationMatrix(wings[sD]->x2, x2_gl);
                arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, wings[sD]->z, x1_gl, x2_gl);
                arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat dz_gldx1 = D1_gl*z_gl;
                arma::mat dz_gldx2 = z_gl*D2_gl.t();
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
                arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
                arma::mat xC = wings[tD]->xC;
                arma::mat yC = wings[tD]->yC;
                arma::mat zC = wings[tD]->zC;
                bw(tD, sD).set_size(wings[tD]->nxy, wings[sD]->nxy, wings[sD]->con);
                for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                    for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                    {
                        size_t k = i+j*wings[tD]->nx;
                        for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                            for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
                                double r3 = pow(norm(r), 3);
                                arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                                                dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                                                dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
                                arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                                                {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                                                {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
                                for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                {
                                    double  t2 = boost::math::chebyshev_t(q, x2_gl(jj));
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                    for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                    {
                                        double  t1 = boost::math::chebyshev_t(p, x1_gl(ii));
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                        arma::vec::fixed<2> dmudxi = {dt1 *  t2, t1 * dt2};
                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                        arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3;
                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k));
                                    }
                                }
                            } 
                    }
            }
    }
    if (sym == Symmetry::y)
    {
        for (size_t sD = 0; sD < wings.size(); sD++)
            for (size_t tD = 0; tD < wings.size(); tD++)
                if (sD != tD)
                {
                    std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                    arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                    {
                        gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                        x1_gl(ii) =-gl_x[ii].x();
                    }
                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                    {
                        gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                        x2_gl(jj) =-gl_y[jj].x();
                    }
                    arma::mat Tx = Lagrange::interpolationMatrix(wings[sD]->x1, x1_gl);
                    arma::mat Ty = Lagrange::interpolationMatrix(wings[sD]->x2, x2_gl);
                    arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, wings[sD]->z, x1_gl, x2_gl);
                    arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
                    arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
                    arma::mat dz_gldx1 = D1_gl*z_gl;
                    arma::mat dz_gldx2 = z_gl*D2_gl.t();
                    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                    y_gl     =-y_gl;
                    dy_gldx1 =-dy_gldx1;
                    dy_gldx2 =-dy_gldx2;
                    arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
                    arma::cube e_c_gl = MetricCo(J_gl);
                    arma::cube ec_gl  = MetricContra(e_c_gl);
                    arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
                    arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
                    arma::mat xC = wings[tD]->xC;
                    arma::mat yC = wings[tD]->yC;
                    arma::mat zC = wings[tD]->zC;
                    for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                        for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                        {
                            size_t k = i+j*wings[tD]->nx;
                            for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                                {
                                    arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
                                    double r3 = pow(norm(r), 3);
                                    arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                                                    dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                                                    dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
                                    arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                                                    {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                                                    {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
                                    for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    {
                                        double  t2 = boost::math::chebyshev_t(q, x2_gl(jj));
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                        for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        {
                                            double  t1 = boost::math::chebyshev_t(p, x1_gl(ii));
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                            arma::vec::fixed<2> dmudxi = {dt1 *  t2, t1 * dt2};
                                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                            arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                                            arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3;
                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k));
                                        }
                                    }
                                } 
                        }
                }
    }
    solve();
}

void Aerodynamics::solve()
{
    bool converged = false;
    int count = 1;
    arma::field<arma::vec> muTarget(interfaces.size()), muSource(interfaces.size());
    arma::field<arma::mat> b0(wings.size());
    for (size_t w = 0; w < wings.size(); w++)
        b0(w) = wings[w]->b;
    for (auto& interface:interfaces)
    {
        Direction targetDirection = static_cast<Direction>(interface.targetCurve);
        Wing *wingTarget = wings[interface.targetDomain];
        switch (interface.sourceCurve)
        {
            case 0: // South
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        wingTarget->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
        }
        Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
        Wing *wingSource = wings[interface.sourceDomain];
        switch (interface.targetCurve)
        {
            case 0: // South
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        wingSource->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
        }
    }
    // Calculation for each wing surface
    switch (analysis)
    {
        case Analysis::linear:
            #pragma omp parallel for
            for (auto& wing:wings)
                wing->linearSolve();
            break;
        case Analysis::nonlinear:
            #pragma omp parallel for
            for (auto& wing:wings)
                wing->nonlinearSolve();
            break;
    }
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (auto& interface:interfaces)
        {
            Wing *wingSource = wings[interface.sourceDomain];
            Wing *wingTarget = wings[interface.targetDomain];
            size_t nx = wingSource->nx;
            size_t ny = wingSource->ny;
            arma::mat mu_hat = wingSource->mu_hat;
            arma::mat  T1 = wingSource->T1;
            arma::mat  T2 = wingSource->T2;
            arma::mat dT1 = wingSource->dT1;
            arma::mat dT2 = wingSource->dT2;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double  t2 = pow(-1, q);
                        double dt2 = pow(-1, q+1) * pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingSource->h_2s1_south;
                    arma::vec h_2s2 = wingSource->h_2s2_south;
                    arma::vec sourceMU = interface.lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            sourceMU = reverse(sourceMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i) = sourceMU(i);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[0]->curveType == CurveType::Boundary  && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1) = sourceMU(nx-1);
                            break;
                        case 1: // East
                            sourceMU = reverse(sourceMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(wingTarget->nx-1+j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(wingTarget->nx-1) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(wingTarget->nx-1+(ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i+(wingTarget->ny-1)*nx) = sourceMU(i);
                            if (wingTarget->chi[2]->curveType == CurveType::Boundary  && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((wingTarget->ny-1)*nx) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1+(wingTarget->ny-1)*nx) = sourceMU(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double dt1 = pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *      dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingSource->h_1s1_east;
                    arma::rowvec h_1s2 = wingSource->h_1s2_east;
                    arma::vec sourceMU = (interface.lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            sourceMU = reverse(sourceMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i) = sourceMU(i);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[0]->curveType == CurveType::Boundary  && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1) = sourceMU(nx-1);
                            break;
                        case 1: // East
                            sourceMU = reverse(sourceMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(wingTarget->nx-1+j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(wingTarget->nx-1) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(wingTarget->nx-1+(ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i+(wingTarget->ny-1)*nx) = sourceMU(i);
                            if (wingTarget->chi[2]->curveType == CurveType::Boundary  && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((wingTarget->ny-1)*nx) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1+(wingTarget->ny-1)*nx) = sourceMU(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double dt2 = pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingSource->h_2s1_north;
                    arma::vec h_2s2 = wingSource->h_2s2_north;
                    arma::vec sourceMU = interface.lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i) = sourceMU(i);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[0]->curveType == CurveType::Boundary  && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1) = sourceMU(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(wingTarget->nx-1+j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(wingTarget->nx-1) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(wingTarget->nx-1+(ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                        case 2: // North
                            sourceMU = reverse(sourceMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i+(wingTarget->ny-1)*nx) = sourceMU(i);
                            if (wingTarget->chi[2]->curveType == CurveType::Boundary  && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((wingTarget->ny-1)*nx) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1+(wingTarget->ny-1)*nx) = sourceMU(nx-1);
                            break;
                        case 3: // West
                            sourceMU = reverse(sourceMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double  t1 = pow(-1, p);
                        double dt1 = pow(-1, p+1) * pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *  t1 * dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingSource->h_1s1_west;
                    arma::rowvec h_1s2 = wingSource->h_1s2_west;
                    arma::vec sourceMU = (interface.lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i) = sourceMU(i);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[0]->curveType == CurveType::Boundary  && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1) = sourceMU(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(wingTarget->nx-1+j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(wingTarget->nx-1) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(wingTarget->nx-1+(ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                        case 2: // North
                            sourceMU = reverse(sourceMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.targetDomain)(i+(wingTarget->ny-1)*nx) = sourceMU(i);
                            if (wingTarget->chi[2]->curveType == CurveType::Boundary  && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((wingTarget->ny-1)*nx) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[1]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)(nx-1+(wingTarget->ny-1)*nx) = sourceMU(nx-1);
                            break;
                        case 3: // West
                            sourceMU = reverse(sourceMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.targetDomain)(j*wingTarget->nx) = sourceMU(j);
                            if (wingTarget->chi[0]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.targetDomain)(0) = sourceMU(0);
                            if (wingTarget->chi[2]->curveType == CurveType::Interface && wingTarget->chi[3]->curveType == CurveType::Interface)
                                b0(interface.targetDomain)((ny-1)*wingTarget->nx) = sourceMU(ny-1);
                            break;
                    }
                    break;
                }
            }
            nx     = wingTarget->nx;
            ny     = wingTarget->ny;
            mu_hat = wingTarget->mu_hat;
             T1    = wingTarget->T1;
             T2    = wingTarget->T2;
            dT1    = wingTarget->dT1;
            dT2    = wingTarget->dT2;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double t2  = pow(-1, q);
                        double dt2 = pow(-1, q+1) * pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingTarget->h_2s1_south;
                    arma::vec h_2s2 = wingTarget->h_2s2_south;
                    arma::vec targetMU = interface.lambdaTarget*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            targetMU = reverse(targetMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i) = targetMU(i);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[0]->curveType == CurveType::Boundary  && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1) = targetMU(nx-1);
                            break;
                        case 1: // East
                            targetMU = reverse(targetMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(wingSource->nx-1+j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(wingSource->nx-1) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(wingSource->nx-1+(ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i+(wingSource->ny-1)*nx) = targetMU(i);
                            if (wingSource->chi[2]->curveType == CurveType::Boundary  && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((wingSource->ny-1)*nx) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1+(wingSource->ny-1)*nx) = targetMU(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double dt1 = pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *      dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_east;
                    arma::vec targetMU = (interface.lambdaTarget*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            targetMU = reverse(targetMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i) = targetMU(i);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[0]->curveType == CurveType::Boundary  && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1) = targetMU(nx-1);
                            break;
                        case 1: // East
                            targetMU = reverse(targetMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(wingSource->nx-1+j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(wingSource->nx-1) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(wingSource->nx-1+(ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i+(wingSource->ny-1)*nx) = targetMU(i);
                            if (wingSource->chi[2]->curveType == CurveType::Boundary  && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((wingSource->ny-1)*nx) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1+(wingSource->ny-1)*nx) = targetMU(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double dt2 = pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wings[interface.targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interface.targetDomain]->h_2s2_north;
                    arma::vec targetMU = interface.lambdaTarget*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i) = targetMU(i);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[0]->curveType == CurveType::Boundary  && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1) = targetMU(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(wingSource->nx-1+j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(wingSource->nx-1) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(wingSource->nx-1+(ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                        case 2: // North
                            targetMU = reverse(targetMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i+(wingSource->ny-1)*nx) = targetMU(i);
                            if (wingSource->chi[2]->curveType == CurveType::Boundary  && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((wingSource->ny-1)*nx) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1+(wingSource->ny-1)*nx) = targetMU(nx-1);
                            break;
                        case 3: // West
                            targetMU = reverse(targetMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double  t1 = pow(-1, p);
                        double dt1 = pow(-1, p+1) * pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {

                                MU(j)    += mu_hat(p+q*nx, 0) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *  t1 * dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_west;
                    arma::vec targetMU = (interface.lambdaTarget*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i) = targetMU(i);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[0]->curveType == CurveType::Boundary  && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1) = targetMU(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(wingSource->nx-1+j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(wingSource->nx-1+(ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                        case 2: // North
                            targetMU = reverse(targetMU);
                            for (size_t i = 1; i < nx-1; i++)
                                b0(interface.sourceDomain)(i+(wingSource->ny-1)*nx) = targetMU(i);
                            if (wingSource->chi[2]->curveType == CurveType::Boundary  && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((wingSource->ny-1)*nx) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[1]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)(nx-1+(wingSource->ny-1)*nx) = targetMU(nx-1);
                            break;
                        case 3: // West
                            targetMU = reverse(targetMU);
                            for (size_t j = 1; j < ny-1; j++)
                                b0(interface.sourceDomain)(j*wingSource->nx) = targetMU(j);
                            if (wingSource->chi[0]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Boundary)
                                b0(interface.sourceDomain)(0) = targetMU(0);
                            if (wingSource->chi[2]->curveType == CurveType::Interface && wingSource->chi[3]->curveType == CurveType::Interface)
                                b0(interface.sourceDomain)((ny-1)*wingSource->nx) = targetMU(ny-1);
                            break;
                    }
                    break;
                }
            }
        }
        // Calculate the influence of the wing surfaces on each other
        for (size_t tD = 0; tD < wings.size(); tD++)
            for (size_t sD = 0; sD < wings.size(); sD++)
                if (sD != tD)
                    wings[tD]->b = b0(tD) + bw(tD, sD).slice(0)*wings[sD]->mu_hat;
        // Calculation for each wing surface
        switch (analysis)
        {
            case Analysis::linear:
                #pragma omp parallel for
                for (auto& wing:wings)
                    wing->linearEval();
                break;
            case Analysis::nonlinear:
                #pragma omp parallel for
                for (auto& wing:wings)
                    wing->nonlinearEval();
                break;
        }

        for (size_t k = 0; k < interfaces.size(); k++)
        {
            Interface interface = interfaces[k];
            Wing *wingSource = wings[interface.sourceDomain];
            size_t nx = wingSource->nx;
            size_t ny = wingSource->ny;
            arma::mat mu_hat = wingSource->mu_hat;
            arma::mat  T1 = wingSource->T1;
            arma::mat  T2 = wingSource->T2;
            arma::mat dT1 = wingSource->dT1;
            arma::mat dT2 = wingSource->dT2;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double  t2 = pow(-1, q);
                        double dt2 = pow(-1, q+1) * pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingSource->h_2s1_south;
                    arma::vec h_2s2 = wingSource->h_2s2_south;
                    muSource(k) = interface.lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double dt1 = pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {                                
                                MU(j)    += mu_hat(p+q*nx, 0) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *      dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingSource->h_1s1_east;
                    arma::rowvec h_1s2 = wingSource->h_1s2_east;
                    muSource(k) = (interface.lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double dt2 = pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingSource->h_2s1_north;
                    arma::vec h_2s2 = wingSource->h_2s2_north;
                    muSource(k) = interface.lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction   
                    {
                        double  t1 = pow(-1, p);
                        double dt1 = pow(-1, p+1) * pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *  t1 * dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingSource->h_1s1_west;
                    arma::rowvec h_1s2 = wingSource->h_1s2_west;
                    muSource(k) = (interface.lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    break;
                }
            }
            Wing *wingTarget = wings[interface.targetDomain];
            nx = wingTarget->nx;
            ny = wingTarget->ny;
            mu_hat = wingTarget->mu_hat;
             T1 = wingTarget->T1;
             T2 = wingTarget->T2;
            dT1 = wingTarget->dT1;
            dT2 = wingTarget->dT2;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double  t2 = pow(-1, q);
                        double dt2 = pow(-1, q+1) * pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p) * t2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p) * t2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingTarget->h_2s1_south;
                    arma::vec h_2s2 = wingTarget->h_2s2_south;
                    muTarget(k) = interface.lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double dt1 = pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *      dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingTarget->h_1s1_east;
                    arma::rowvec h_1s2 = wingTarget->h_1s2_east;
                    muTarget(k) = (interface.lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                    {
                        double dt2 = pow(q, 2);
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                            {
                                MU(i)    += mu_hat(p+q*nx, 0) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx, 0) *  T1(i, p) * dt2;
                            }
                    }
                    arma::vec h_2s1 = wingTarget->h_2s1_north;
                    arma::vec h_2s2 = wingTarget->h_2s2_north;
                    muTarget(k) = interface.lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                    {
                        double  t1 = pow(-1, p);
                        double dt1 = pow(-1, p+1) * pow(p, 2);
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                            {
                                MU(j)    += mu_hat(p+q*nx, 0) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx, 0) *  t1 * dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingTarget->h_1s1_west;
                    arma::rowvec h_1s2 = wingTarget->h_1s2_west;
                    muTarget(k) = (interface.lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
            }
        }
        // Convergence check
        converged = true;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            printf("Interface %lu: ", k+1);
            double res = fabs(1 - arma::norm(muSource(k).subvec(1, muSource(k).size()-2))/arma::norm(muTarget(k).subvec(1, muTarget(k).size()-2)));
            printf("Residual %4.2e\n", res);
            if (res > residualTarget)
                converged = false;
        }
        count++;
        if (count > iterations)
            converged = true;
    } while (converged == false);
    #pragma omp parallel for
    for (auto& wing:wings)
        wing->postprocessing();
}

arma::vec Aerodynamics::get_lift()
{
    arma::vec lift(wings[0]->con, arma::fill::zeros);
    for (auto wing:wings)
        lift += wing->lift;
    return lift;
}

arma::vec Aerodynamics::get_moment()
{
    arma::vec moment(wings[0]->con, arma::fill::zeros);
    for (auto wing:wings)
        moment += wing->moment;
    return moment;
}
double Aerodynamics::get_area()
{
    double area = 0;
    for (auto wing:wings)
        area += wing->area;
    return area;
}

void Aerodynamics::output(const std::string filename)
{
    for (size_t k = 0; k < wings.size(); k++)
        wings[k]->output(filename+"_"+std::to_string(k));
}

template void Aerodynamics::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const int);
template void Aerodynamics::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const size_t);
template void Aerodynamics::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double);
template void Aerodynamics::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const arma::vec);
template void Aerodynamics::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const int);
template void Aerodynamics::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const size_t);
template void Aerodynamics::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const double);
template void Aerodynamics::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const arma::vec);