#include "Aerodynamics.hpp"

void Aerodynamics::checkMesh()
{
    for (size_t i = 0; i < wings.size(); i++)
    {
        printf("Checking mesh of wing number %lu\n", i);
        wings[i]->checkMesh();
    }
}

void Aerodynamics::output(const std::string filename)
{
    for (size_t k = 0; k < wings.size(); k++)
    {
        std::string file_k = filename;
        file_k.append("_");
        file_k.append(std::to_string(k));
        wings[k]->output(file_k);
    }
}

void Aerodynamics::linear()
{
    analysis = Analysis::linear;
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

    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (auto& interface:interfaces)
        {
            Direction targetDirection = static_cast<Direction>(interface.targetCurve);
            Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
            size_t nx = wings[interface.sourceDomain]->nx;
            size_t ny = wings[interface.sourceDomain]->ny;
            arma::vec x1 = wings[interface.sourceDomain]->x1;
            arma::vec x2 = wings[interface.sourceDomain]->x2;
            arma::mat MU     = wings[interface.sourceDomain]->mu;
            arma::mat mu_hat = wings[interface.sourceDomain]->mu_hat;
            arma::mat dMUd1(nx, ny, arma::fill::zeros);
            arma::mat dMUd2(nx, ny, arma::fill::zeros);
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
                            dMUd1(i, j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                            dMUd2(i, j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                        }
                    }
            double gamma;
            arma::vec sourceMU;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec h_2s1 = wings[interface.sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interface.sourceDomain]->h_2s2_south;
                    gamma = gamma0*mean(h_2s2);
                    sourceMU = gamma*MU.col(0) + h_2s1%dMUd1.col(0) + h_2s2%dMUd2.col(0);
                    break;
                }
                case 1: // East
                {
                    arma::rowvec h_1s1 = wings[interface.sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interface.sourceDomain]->h_1s2_east;
                    gamma = gamma0*mean(h_1s1);
                    sourceMU = (gamma*MU.row(nx-1) - (h_1s1%dMUd1.row(nx-1) + h_1s2%dMUd2.row(nx-1))).t();
                    break;
                }
                case 2: // North
                {
                    arma::vec h_2s1 = wings[interface.sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interface.sourceDomain]->h_2s2_north;
                    gamma = gamma0*mean(h_2s2);
                    sourceMU = gamma*MU.col(ny-1) - (h_2s1%dMUd1.col(ny-1) + h_2s2%dMUd2.col(ny-1));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = wings[interface.sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interface.sourceDomain]->h_1s2_west;
                    gamma = gamma0*mean(h_1s1);
                    sourceMU = (gamma*MU.row(0) + h_1s1%dMUd1.row(0) + h_1s2%dMUd2.row(0)).t();
                    break;
                }
            }
            switch (interface.targetCurve)
            {
                case 0: case 3: // South and West
                    wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, gamma,-1, sourceMU);
                    break;
                case 1: case 2: // East and North
                    wings[interface.targetDomain]->boundary(targetDirection, BC::Robin, gamma, 1, sourceMU);
                    break;
            }
            nx = wings[interface.targetDomain]->nx;
            ny = wings[interface.targetDomain]->ny;
            x1 = wings[interface.targetDomain]->x1;
            x2 = wings[interface.targetDomain]->x2;
            MU = wings[interface.targetDomain]->mu;
            mu_hat = wings[interface.targetDomain]->mu_hat;
            dMUd1.zeros(nx, ny);
            dMUd2.zeros(nx, ny);
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
                            dMUd1(i, j) += mu_hat(p+q*nx, 0) * dT1dx1 * T2;
                            dMUd2(i, j) += mu_hat(p+q*nx, 0) * T1 * dT2dx2;
                        }
                    }
            arma::vec targetMU;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec h_2s1 = wings[interface.targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = wings[interface.targetDomain]->h_2s2_south;
                    gamma = gamma0*mean(h_2s2);
                    targetMU = gamma*MU.col(0) + h_2s1%dMUd1.col(0) + h_2s2%dMUd2.col(0);
                    break;
                }
                case 1: // East
                {
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_east;
                    gamma = gamma0*mean(h_1s1);
                    targetMU = (gamma*MU.row(nx-1) - (h_1s1%dMUd1.row(nx-1) + h_1s2%dMUd2.row(nx-1))).t();
                    break;
                }
                case 2: // North
                {
                    arma::vec h_2s1 = wings[interface.targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = wings[interface.targetDomain]->h_2s2_north;
                    gamma = gamma0*mean(h_2s2);
                    targetMU = gamma*MU.col(ny-1) - (h_2s1%dMUd1.col(ny-1) + h_2s2%dMUd2.col(ny-1));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = wings[interface.targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = wings[interface.targetDomain]->h_1s2_west;
                    gamma = gamma0*mean(h_1s1);
                    targetMU = (gamma*MU.row(0) + h_1s1%dMUd1.row(0) + h_1s2%dMUd2.row(0)).t();
                    break;
                }
            }
            switch (interface.sourceCurve)
            {
                case 0: case 3: // South and West
                    wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, gamma,-1, targetMU);
                    break;
                case 1: case 2: // East and North
                    wings[interface.sourceDomain]->boundary(sourceDirection, BC::Robin, gamma, 1, targetMU);
                    break;
            }
        }
        switch (analysis)
        {
            case Analysis::linear:
                #pragma omp parallel for
                for (auto& wing:wings)
                    wing->linear();
                break;
            case Analysis::nonlinear:
                #pragma omp parallel for
                for (auto& wing:wings)
                    wing->nonlinear();
                break;
        }
        count++;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            switch (interfaces[k].targetCurve)
            {
                case 0: // South
                {
                    muTarget(k) = wings[interfaces[k].targetDomain]->mu.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = wings[interfaces[k].targetDomain]->nx;
                    muTarget(k) = wings[interfaces[k].targetDomain]->mu.row(nx-1).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = wings[interfaces[k].targetDomain]->ny;
                    muTarget(k) = wings[interfaces[k].targetDomain]->mu.col(ny-1);
                    break;
                }
                case 3: // West
                {
                    muTarget(k) = wings[interfaces[k].targetDomain]->mu.row(0).t();
                    break;
                }
            }
            switch (interfaces[k].sourceCurve)
            {
                case 0: // South
                {
                    muSource(k) = wings[interfaces[k].sourceDomain]->mu.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = wings[interfaces[k].sourceDomain]->nx;
                    muSource(k) = wings[interfaces[k].sourceDomain]->mu.row(nx-1).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = wings[interfaces[k].sourceDomain]->ny;
                    muSource(k) = wings[interfaces[k].sourceDomain]->mu.col(ny-1);
                    break;
                }
                case 3: // West
                {
                    muSource(k) = wings[interfaces[k].sourceDomain]->mu.row(0).t();
                    break;
                }
            }
        }
        converged = true;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            printf("Interface %lu\n", k+1);
            double res = arma::norm(muTarget(k) - muSource(k))/arma::norm(muTarget(k) + muSource(k));
            printf("Residual %4.2e\n", res);
            if (res > residualTarget)
                converged = false;
        }
        if (count > iterations)
            converged = true;
    } while (converged == false);
}