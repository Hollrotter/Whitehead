#include "Aerodynamics.hpp"

Aerodynamics::Aerodynamics(std::vector<Wing*> _wings) : wings(_wings)
{
    for (size_t sD = 0; sD < wings.size()-1; sD++)
        for (size_t sC = 0; sC < 4; sC++)
            for (size_t tD = sD+1; tD < wings.size(); tD++)
                for (size_t tC = 0; tC < 4 ; tC++)
                    if (wings[sD]->chi[sC] == wings[tD]->chi[tC])
                    {
                        interfaces.push_back(Interface(sD, tD, sC, tC));
                        wings[sD]->chi[sC]->curveType = CurveType::Interface;
                    }
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
                                    double T2     = boost::math::chebyshev_t(q, x2_gl(jj));
                                    double dT2dx2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                    for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction of source
                                    {
                                        double T1     = boost::math::chebyshev_t(p, x1_gl(ii));
                                        double dT1dx1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                        double dmudx1 = dT1dx1*T2;
                                        double dmudx2 = T1*dT2dx2;
                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= gl_x[ii].weight * gl_y[jj].weight / r3
                                            *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                            - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
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
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (auto& interface:interfaces)
        {
            Direction targetDirection = static_cast<Direction>(interface.targetCurve);
            Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
            size_t nx = wings[interface.sourceDomain]->nx;
            size_t ny = wings[interface.sourceDomain]->ny;
            arma::vec xi_1 = wings[interface.sourceDomain]->xi_1;
            arma::vec xi_2 = wings[interface.sourceDomain]->xi_2;
            arma::mat mu_hat = wings[interface.sourceDomain]->mu_hat;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));  
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double T2     = pow(-1, q);
                                double dT2dx2 = pow(-1, q+1) * pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interface.sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interface.sourceDomain]->h_2s2_south;
                    arma::vec sourceMU = interface.lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceMU)));
                            break;
                        case 1: // East
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceMU)));
                            break;
                        case 2: // North
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, sourceMU);
                            break;
                        case 3: // West
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, sourceMU);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double dT1dx1 = pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interface.sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interface.sourceDomain]->h_1s2_east;
                    arma::vec sourceMU = (interface.lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceMU)));
                            break;
                        case 1: // East
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceMU)));
                            break;
                        case 2: // North
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, sourceMU);
                            break;
                        case 3: // West
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, sourceMU);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double T1     = boost::math::chebyshev_t(p, xi_1(i));
                                double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                                for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                {
                                    double dT2dx2 = pow(q, 2);
                                    MU(i)    += mu_hat(p+q*nx, 0) * T1;
                                    dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1;
                                    dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                                }
                            }
                    arma::vec h_2s1 = wings[interface.sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interface.sourceDomain]->h_2s2_north;
                    arma::vec sourceMU = interface.lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, sourceMU);
                            break;
                        case 1: // East
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, sourceMU);
                            break;
                        case 2: // North
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceMU)));
                            break;
                        case 3: // West
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceMU)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction        
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double T1     = pow(-1, p);
                                double dT1dx1 = pow(-1, p+1) * pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interface.sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interface.sourceDomain]->h_1s2_west;
                    arma::vec sourceMU = (interface.lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, sourceMU);
                            break;
                        case 1: // East
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, sourceMU);
                            break;
                        case 2: // North
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceMU)));
                            break;
                        case 3: // West
                            wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceMU)));
                            break;
                    }
                    break;
                }
            }
            nx     = wings[interface.targetDomain]->nx;
            ny     = wings[interface.targetDomain]->ny;
            xi_1   = wings[interface.targetDomain]->xi_1;
            xi_2   = wings[interface.targetDomain]->xi_2;
            mu_hat = wings[interface.targetDomain]->mu_hat;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double T2     = pow(-1, q);
                                double dT2dx2 = pow(-1, q+1) * pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interface.targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interface.targetDomain]->h_2s2_south;
                    arma::vec targetMU = interface.lambdaTarget*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetMU)));
                            break;
                        case 1: // East
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetMU)));
                            break;
                        case 2: // North
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetMU);
                            break;
                        case 3: // West
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetMU);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double dT1dx1 = pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_east;
                    arma::vec targetMU = (interface.lambdaTarget*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetMU)));
                            break;
                        case 1: // East
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetMU)));
                            break;
                        case 2: // North
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetMU);
                            break;
                        case 3: // West
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetMU);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction  
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double dT2dx2 = pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interface.targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interface.targetDomain]->h_2s2_north;
                    arma::vec targetMU = interface.lambdaTarget*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetMU);
                            break;
                        case 1: // East
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetMU);
                            break;
                        case 2: // North
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetMU)));
                            break;
                        case 3: // West
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetMU)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double T1     = pow(-1, p);
                                double dT1dx1 = pow(-1, p+1) * pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_west;
                    arma::vec targetMU = (interface.lambdaTarget*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetMU);
                            break;
                        case 1: // East
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetMU);
                            break;
                        case 2: // North
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetMU)));
                            break;
                        case 3: // West
                            wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetMU)));
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
                // #pragma omp parallel for
                for (auto& wing:wings)
                    wing->linear();
                break;
            case Analysis::nonlinear:
                #pragma omp parallel for
                for (auto& wing:wings)
                    wing->nonlinear();
                break;
        }

        for (size_t k = 0; k < interfaces.size(); k++)
        {
            switch (interfaces[k].sourceCurve)
            {
                case 0: // South
                {
                    size_t nx = wings[interfaces[k].sourceDomain]->nx;
                    size_t ny = wings[interfaces[k].sourceDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].sourceDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].sourceDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].sourceDomain]->mu_hat;
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction 
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double T2     = pow(-1, q);
                                double dT2dx2 = pow(-1, q+1) * pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interfaces[k].sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interfaces[k].sourceDomain]->h_2s2_south;
                    muSource(k) = interfaces[k].lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    break;
                }
                case 1: // East
                {
                    size_t nx = wings[interfaces[k].sourceDomain]->nx;
                    size_t ny = wings[interfaces[k].sourceDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].sourceDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].sourceDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].sourceDomain]->mu_hat;
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double dT1dx1 = pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interfaces[k].sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interfaces[k].sourceDomain]->h_1s2_east;
                    muSource(k) = (interfaces[k].lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    break;
                }
                case 2: // North
                {
                    size_t nx = wings[interfaces[k].sourceDomain]->nx;
                    size_t ny = wings[interfaces[k].sourceDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].sourceDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].sourceDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].sourceDomain]->mu_hat;
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double dT2dx2 = pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interfaces[k].sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interfaces[k].sourceDomain]->h_2s2_north;
                    muSource(k) = interfaces[k].lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    break;
                }
                case 3: // West
                {
                    size_t nx = wings[interfaces[k].sourceDomain]->nx;
                    size_t ny = wings[interfaces[k].sourceDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].sourceDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].sourceDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].sourceDomain]->mu_hat;
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction      
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double T1     = pow(-1, p);
                                double dT1dx1 = pow(-1, p+1) * pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interfaces[k].sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interfaces[k].sourceDomain]->h_1s2_west;
                    muSource(k) = (interfaces[k].lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    break;
                }
            }
            switch (interfaces[k].targetCurve)
            {
                case 0: // South
                {
                    size_t nx = wings[interfaces[k].targetDomain]->nx;
                    size_t ny = wings[interfaces[k].targetDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].targetDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].targetDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].targetDomain]->mu_hat;
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double T2     = pow(-1, q);
                                double dT2dx2 = pow(-1, q+1) * pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interfaces[k].targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interfaces[k].targetDomain]->h_2s2_south;
                    muTarget(k) = interfaces[k].lambdaSource*MU - (h_2s1%dMUd1 + h_2s2%dMUd2);
                    if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 1: // East
                {
                    size_t nx = wings[interfaces[k].targetDomain]->nx;
                    size_t ny = wings[interfaces[k].targetDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].targetDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].targetDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].targetDomain]->mu_hat;
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction 
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double dT1dx1 = pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interfaces[k].targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interfaces[k].targetDomain]->h_1s2_east;
                    muTarget(k) = (interfaces[k].lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 2: // North
                {
                    size_t nx = wings[interfaces[k].targetDomain]->nx;
                    size_t ny = wings[interfaces[k].targetDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].targetDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].targetDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].targetDomain]->mu_hat;
                    arma::vec    MU(nx, arma::fill::zeros);
                    arma::vec dMUd1(nx, arma::fill::zeros);
                    arma::vec dMUd2(nx, arma::fill::zeros);
                    for (size_t i = 0; i < nx; i++) // Loop over Collocation Points in 1-direction
                        for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                        {
                            double T1     = boost::math::chebyshev_t(p, xi_1(i));
                            double dT1dx1 = boost::math::chebyshev_t_prime(p, xi_1(i));
                            for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                            {
                                double dT2dx2 = pow(q, 2);
                                MU(i)    += mu_hat(p+q*nx, 0) * T1;
                                dMUd1(i) += mu_hat(p+q*nx, 0) * dT1dx1;
                                dMUd2(i) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::vec h_2s1 = wings[interfaces[k].targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interfaces[k].targetDomain]->h_2s2_north;
                    muTarget(k) = interfaces[k].lambdaSource*MU + h_2s1%dMUd1 + h_2s2%dMUd2;
                    if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
                        muTarget(k) = reverse(muTarget(k));
                    break;
                }
                case 3: // West
                {
                    size_t nx = wings[interfaces[k].targetDomain]->nx;
                    size_t ny = wings[interfaces[k].targetDomain]->ny;
                    arma::vec xi_1 = wings[interfaces[k].targetDomain]->xi_1;
                    arma::vec xi_2 = wings[interfaces[k].targetDomain]->xi_2;
                    arma::mat mu_hat = wings[interfaces[k].targetDomain]->mu_hat;
                    arma::rowvec    MU(ny, arma::fill::zeros);
                    arma::rowvec dMUd1(ny, arma::fill::zeros);
                    arma::rowvec dMUd2(ny, arma::fill::zeros);
                    for (size_t j = 0; j < ny; j++) // Loop over Collocation Points in 2-direction
                        for (size_t q = 0; q < ny; q++) // Loop over Chebyshev Polynomial 2-direction
                        {
                            double T2     = boost::math::chebyshev_t(q, xi_2(j));
                            double dT2dx2 = boost::math::chebyshev_t_prime(q, xi_2(j));
                            for (size_t p = 0; p < nx; p++) // Loop over Chebyshev Polynomial 1-direction
                            {
                                double T1     = pow(-1, p);
                                double dT1dx1 = pow(-1, p+1) * pow(p, 2);
                                MU(j)    += mu_hat(p+q*nx, 0) * T1 * T2;
                                dMUd1(j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                                dMUd2(j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                            }
                        }
                    arma::rowvec h_1s1 = wings[interfaces[k].targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interfaces[k].targetDomain]->h_1s2_west;
                    muTarget(k) = (interfaces[k].lambdaSource*MU - (h_1s1%dMUd1 + h_1s2%dMUd2)).t();
                    if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
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