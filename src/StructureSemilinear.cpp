#include "Structure.hpp"

void Structure::semilinear()
{
    analysis = Analysis::semilinear;
    bool converged = false;
    size_t count = 1;
    arma::field<arma::vec> Ztarget(interfaces.size()), Zsource(interfaces.size());

    arma::vec iter_old(membranes.size());
    for (size_t k = 0; k < membranes.size(); k++)
    {
        iter_old(k) = membranes[k]->iter;
        membranes[k]->iter = 1;
    }
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << '\n';
        for (Interface& interface:interfaces)
        {
            Direction targetDirection = static_cast<Direction>(interface.targetCurve);
            Membrane *membraneSource = membranes[interface.sourceDomain];
            Membrane *membraneTarget = membranes[interface.targetDomain];
            arma::mat D1 = membraneSource->D1;
            arma::mat D2 = membraneSource->D2;
            arma::mat Z  = membraneSource->z;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membraneSource->h_2s1_south;
                    arma::vec h_2s2 = membraneSource->h_2s2_south;
                    arma::vec sourceZ = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                        case 1: // East
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 2: // North
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 3: // West
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneSource->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                    arma::vec sourceZ = (interface.lambdaSource*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                        case 1: // East
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 2: // North
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 3: // West
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneSource->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneSource->h_2s1_north;
                    arma::vec h_2s2 = membraneSource->h_2s2_north;
                    arma::vec sourceZ = interface.lambdaSource*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                        case 1: // East
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 2: // North
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 3: // West
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membraneSource->h_1s1_west;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_west;
                    arma::vec sourceZ = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                        case 1: // East
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 2: // North
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 3: // West
                            membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                    }
                    break;
                }
            }
            Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
            D1 = membraneTarget->D1;
            D2 = membraneTarget->D2;
            Z  = membraneTarget->z;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membraneTarget->h_2s1_south;
                    arma::vec h_2s2 = membraneTarget->h_2s2_south;
                    arma::vec targetZ = interface.lambdaTarget*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                        case 1: // East
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 2: // North
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 3: // West
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneTarget->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_east;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_east;
                    arma::vec targetZ = (interface.lambdaTarget*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                        case 1: // East
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 2: // North
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 3: // West
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneTarget->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneTarget->h_2s1_north;
                    arma::vec h_2s2 = membraneTarget->h_2s2_north;
                    arma::vec targetZ = interface.lambdaTarget*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                        case 1: // East
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 2: // North
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 3: // West
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_west;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_west;
                    arma::vec targetZ = (interface.lambdaTarget*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                        case 1: // East
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 2: // North
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 3: // West
                            membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                    }
                    break;
                }
            }
        }
        #pragma omp parallel for
        for (auto& membrane:membranes)
            membrane->semilinear();
        
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            Interface interface = interfaces[k];
            const Membrane *membraneSource = membranes[interface.sourceDomain];
            arma::mat D1 = membraneSource->D1;
            arma::mat D2 = membraneSource->D2;
            arma::mat Z  = membraneSource->z;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membraneSource->h_2s1_south;
                    arma::vec h_2s2 = membraneSource->h_2s2_south;
                    Zsource(k) = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneSource->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                    Zsource(k) = (interface.lambdaSource*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneSource->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneSource->h_2s1_north;
                    arma::vec h_2s2 = membraneSource->h_2s2_north;
                    Zsource(k) = interface.lambdaSource*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membraneSource->h_1s1_west;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_west;
                    Zsource(k) = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    break;
                }
            }
            const Membrane *membraneTarget = membranes[interface.targetDomain];
            D1 = membraneTarget->D1;
            D2 = membraneTarget->D2;
            Z  = membraneTarget->z;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membraneTarget->h_2s1_south;
                    arma::vec h_2s2 = membraneTarget->h_2s2_south;
                    Ztarget(k) = interface.lambdaSource*Z.col(0) - (h_2s1%dZd1 + h_2s2%dZd2);
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneTarget->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_east;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_east;
                    Ztarget(k) = (interface.lambdaSource*Z.row(nx-1) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneTarget->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneTarget->h_2s1_north;
                    arma::vec h_2s2 = membraneTarget->h_2s2_north;
                    Ztarget(k) = interface.lambdaSource*Z.col(ny-1) + h_2s1%dZd1 + h_2s2%dZd2;
                    if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_west;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_west;
                    Ztarget(k) = (interface.lambdaSource*Z.row(0) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
            }
        }
        converged = true;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            printf("Interface %lu\n", k+1);
            double res = fabs(1 - arma::norm(Zsource(k).subvec(1, Zsource(k).size()-2))/arma::norm(Ztarget(k).subvec(1, Ztarget(k).size()-2)));
            printf("Residual %4.2e\n", res);
            if (res > residualTarget)
                converged = false;
        }
        count++;
        if (count > iterations)
            converged = true;
    } while (converged == false);
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->iter = iter_old(k);
}