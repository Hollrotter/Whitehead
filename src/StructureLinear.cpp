#include "Structure.hpp"

void Structure::linear()
{
    analysis = Analysis::linear;
    bool converged = false;
    int count = 0;
    arma::field<arma::vec> Ztarget(interfaces.size()), Zsource(interfaces.size());
    for (size_t k = 0; k < interfaces.size(); k++)
    {
        size_t nn;
        switch (interfaces[k].sourceCurve)
        {
            case 0: case 2:
                nn = membranes[interfaces[k].sourceDomain]->nx;
                break;
            case 1: case 3:
                nn = membranes[interfaces[k].sourceDomain]->ny;
                break;
        }
    }
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << '\n';
        for (auto& interface:interfaces)
        {
            Direction targetDirection = static_cast<Direction>(interface.targetCurve);
            Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
            arma::mat D1   = membranes[interface.sourceDomain]->D1;
            arma::mat D2   = membranes[interface.sourceDomain]->D2;
            arma::mat Z    = membranes[interface.sourceDomain]->z;
            arma::mat dZd1 = D1*Z;
            arma::mat dZd2 = Z*D2.t();
            double gamma;
            arma::vec sourceZ;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_south;
                    gamma = gamma0*mean(h_2s2);
                    sourceZ = gamma*Z.col(0) + h_2s1%dZd1.col(0) + h_2s2%dZd2.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interface.sourceDomain]->nx;
                    arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_east;
                    gamma = gamma0*mean(h_1s1);
                    sourceZ = (gamma*Z.row(nx-1) - (h_1s1%dZd1.row(nx-1) + h_1s2%dZd2.row(nx-1))).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interface.sourceDomain]->ny;
                    arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_north;
                    gamma = gamma0*mean(h_2s2);
                    sourceZ = gamma*Z.col(ny-1) - (h_2s1%dZd1.col(ny-1) + h_2s2%dZd2.col(ny-1));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_west;
                    gamma = gamma0*mean(h_1s1);
                    sourceZ = (gamma*Z.row(0) + h_1s1%dZd1.row(0) + h_1s2%dZd2.row(0)).t();
                    break;
                }
            }
            switch (interface.targetCurve)
            {
                case 0: case 3: // South and West
                    membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, gamma,-1, sourceZ);
                    break;
                case 1: case 2: // East and North
                    membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, gamma, 1, sourceZ);
                    break;
            }
            D1   = membranes[interface.targetDomain]->D1;
            D2   = membranes[interface.targetDomain]->D2;
            Z    = membranes[interface.targetDomain]->z;
            dZd1 = D1*Z;
            dZd2 = Z*D2.t();
            arma::vec targetZ;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_south;
                    gamma = gamma0*mean(h_2s2);
                    targetZ = gamma*Z.col(0) + h_2s1%dZd1.col(0) + h_2s2%dZd2.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interface.targetDomain]->nx;
                    arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_east;
                    gamma = gamma0*mean(h_1s1);
                    targetZ = (gamma*Z.row(nx-1) - (h_1s1%dZd1.row(nx-1) + h_1s2%dZd2.row(nx-1))).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interface.targetDomain]->ny;
                    arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_north;
                    gamma = gamma0*mean(h_2s2);
                    targetZ = gamma*Z.col(ny-1) - (h_2s1%dZd1.col(ny-1) + h_2s2%dZd2.col(ny-1));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_west;
                    gamma = gamma0*mean(h_1s1);
                    targetZ = (gamma*Z.row(0) + h_1s1%dZd1.row(0) + h_1s2%dZd2.row(0)).t();
                    break;
                }
            }
            switch (interface.sourceCurve)
            {
                case 1: case 2: // East and North
                    membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, gamma, 1, targetZ);
                    break;
                case 0: case 3: // South and West
                    membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, gamma,-1, targetZ);
                    break;
            }
        }
        #pragma omp parallel for
        for (auto& membrane:membranes)
            membrane->linear();
        count++;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            switch (interfaces[k].targetCurve)
            {
                case 0: // South
                {
                    Ztarget(k) = membranes[interfaces[k].targetDomain]->z.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interfaces[k].targetDomain]->nx;
                    Ztarget(k) = membranes[interfaces[k].targetDomain]->z.row(nx-1).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interfaces[k].targetDomain]->ny;
                    Ztarget(k) = membranes[interfaces[k].targetDomain]->z.col(ny-1);
                    break;
                }
                case 3: // West
                {
                    Ztarget(k) = membranes[interfaces[k].targetDomain]->z.row(0).t();
                    break;
                }
            }
            switch (interfaces[k].sourceCurve)
            {
                case 0: // South
                {
                    Zsource(k) = membranes[interfaces[k].sourceDomain]->z.col(0);
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interfaces[k].sourceDomain]->nx;
                    Zsource(k) = membranes[interfaces[k].sourceDomain]->z.row(nx-1).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interfaces[k].sourceDomain]->ny;
                    Zsource(k) = membranes[interfaces[k].sourceDomain]->z.col(ny-1);
                    break;
                }
                case 3: // West
                {
                    Zsource(k) = membranes[interfaces[k].sourceDomain]->z.row(0).t();
                    break;
                }
            }
        }
        converged = true;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            printf("Interface %lu\n", k+1);
            double res = arma::norm(Ztarget(k) - Zsource(k))/arma::norm(Ztarget(k) + Zsource(k));
            printf("Residual %4.2e\n", res);
            if (res > residualTarget)
                converged = false;
        }
        if (count > iterations)
            converged = true;
    } while (converged == false);
}