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
    double l0 = 2;
    for (Interface& interface:_interfaces)
    {
        switch (interface.sourceCurve)
        {
            case 0: // South
                interface.lambdaSource = l0*mean(_wings[interface.sourceDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaSource = l0*mean(_wings[interface.sourceDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaSource = l0*mean(_wings[interface.sourceDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaSource = l0*mean(_wings[interface.sourceDomain]->h_1s1_west);
                break;
        }
        switch (interface.targetCurve)
        {
            case 0: // South
                interface.lambdaTarget = l0*mean(_wings[interface.targetDomain]->h_2s2_south);
                break;
            case 1: // East
                interface.lambdaTarget = l0*mean(_wings[interface.targetDomain]->h_1s1_east);
                break;
            case 2: // North
                interface.lambdaTarget = l0*mean(_wings[interface.targetDomain]->h_2s2_north);
                break;
            case 3: // West
                interface.lambdaTarget = l0*mean(_wings[interface.targetDomain]->h_1s1_west);
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
        std::cout << "Checking mesh of wing number " << i << "\n";
        wings[i]->checkMesh();
    }
}

void Aerodynamics::wake(Wake* w)
{
    for (auto& wing:wings)
        if (std::any_of(wing->chi.begin(), wing->chi.end(), [&](Lagrange::CurveInterpolant* c) {return c == w->chi;}))
        {
            wing->wake(w);
            return;
        }
    std::println("This wake is not part of any of the wing surfaces!");
}

void Aerodynamics::solve()
{
    bool converged = false;
    size_t count = 1;
    arma::field<arma::vec> muTarget(interfaces.size()), muSource(interfaces.size());
    for (const auto& interface:interfaces)
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
        default:
            std::println("Only linear and nonlinear analysis is implemented for Aerodynamics!");
            exit(EXIT_FAILURE);
    }
    arma::field<arma::mat> b0(wings.size());
    for (size_t w = 0; w < wings.size(); w++)
        b0(w) = wings[w]->b;
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (auto& interface:interfaces)
        {
            Wing *wingSource = wings[interface.sourceDomain];
            Wing *wingTarget = wings[interface.targetDomain];
            size_t nx = wingSource->nx;
            size_t ny = wingSource->ny;
            arma::vec mu_hat = wingSource->mu_hat;
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *      dT2(j, q);
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *  t1 * dT2(j, q);
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *      dT2(j, q);
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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

                                MU(j)    += mu_hat(p+q*nx) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *  t1 * dT2(j, q);
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
                    wings[tD]->b = b0(tD) + bw(tD, sD)*wings[sD]->mu_hat;
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
            default:
                std::println("Only linear and nonlinear analysis are implemented for Aerodynamics!");
                exit(EXIT_FAILURE);
        }

        for (size_t k = 0; k < interfaces.size(); k++)
        {
            Interface interface = interfaces[k];
            const Wing *wingSource = wings[interface.sourceDomain];
            size_t nx = wingSource->nx;
            size_t ny = wingSource->ny;
            arma::vec mu_hat = wingSource->mu_hat;
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p) *  t2;
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p) *  t2;
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *      dT2(j, q);
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *  t1 * dT2(j, q);
                            }
                    }
                    arma::rowvec h_1s1 = wingSource->h_1s1_west;
                    arma::rowvec h_1s2 = wingSource->h_1s2_west;
                    muSource(k) = (interface.lambdaSource*MU + h_1s1%dMUd1 + h_1s2%dMUd2).t();
                    break;
                }
            }
            const Wing *wingTarget = wings[interface.targetDomain];
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p) * t2;
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p) * t2;
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *       T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 * T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *      dT2(j, q);
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
                                MU(i)    += mu_hat(p+q*nx) *  T1(i, p);
                                dMUd1(i) += mu_hat(p+q*nx) * dT1(i, p);
                                dMUd2(i) += mu_hat(p+q*nx) *  T1(i, p) * dt2;
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
                                MU(j)    += mu_hat(p+q*nx) *  t1 *  T2(j, q);
                                dMUd1(j) += mu_hat(p+q*nx) * dt1 *  T2(j, q);
                                dMUd2(j) += mu_hat(p+q*nx) *  t1 * dT2(j, q);
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
            std::cout << "Interface " << k+1 << ": ";
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

void Aerodynamics::output(const std::string &filename)
{
    for (size_t k = 0; k < wings.size(); k++)
        wings[k]->output(filename+"_"+std::to_string(k));
}