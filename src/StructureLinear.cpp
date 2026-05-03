#include "Structure.hpp"

void Structure::linear()
{
    analysis = Analysis::linear;
    bool converged = false;
    int count = 1;
    arma::field<arma::vec> Ztarget(interfaces.size()), Zsource(interfaces.size());

    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (Interface& interface:interfaces)
        {
            Direction targetDirection = static_cast<Direction>(interface.targetCurve);
            Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
            arma::mat D1 = membranes[interface.sourceDomain]->D1;
            arma::mat D2 = membranes[interface.sourceDomain]->D2;
            arma::mat Z  = membranes[interface.sourceDomain]->z;
            switch (interface.sourceCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_south;
                    arma::vec sourceZ = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                        case 1: // East
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 2: // North
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 3: // West
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interface.sourceDomain]->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_east;
                    arma::vec sourceZ = (interface.lambdaSource*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                        case 1: // East
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 2: // North
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 3: // West
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interface.sourceDomain]->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_north;
                    arma::vec sourceZ = interface.lambdaSource*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                        case 1: // East
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 2: // North
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 3: // West
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_west;
                    arma::vec sourceZ = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                            break;
                        case 1: // East
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                            break;
                        case 2: // North
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                            break;
                        case 3: // West
                            membranes[interface.targetDomain]->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                            break;
                    }
                    break;
                }
            }
            D1 = membranes[interface.targetDomain]->D1;
            D2 = membranes[interface.targetDomain]->D2;
            Z  = membranes[interface.targetDomain]->z;
            switch (interface.targetCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_south;
                    arma::vec targetZ = interface.lambdaTarget*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                        case 1: // East
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 2: // North
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 3: // West
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interface.targetDomain]->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_east;
                    arma::vec targetZ = (interface.lambdaTarget*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                        case 1: // East
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 2: // North
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 3: // West
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interface.targetDomain]->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_north;
                    arma::vec targetZ = interface.lambdaTarget*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                        case 1: // East
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 2: // North
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 3: // West
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_west;
                    arma::vec targetZ = (interface.lambdaTarget*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                            break;
                        case 1: // East
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                            break;
                        case 2: // North
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                            break;
                        case 3: // West
                            membranes[interface.sourceDomain]->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                            break;
                    }
                    break;
                }
            }
        }
        #pragma omp parallel for
        for (auto& membrane:membranes)
            membrane->linear();

        for (size_t k = 0; k < interfaces.size(); k++)
        {
            arma::mat Z  = membranes[interfaces[k].sourceDomain]->z;
            arma::mat D1 = membranes[interfaces[k].sourceDomain]->D1;
            arma::mat D2 = membranes[interfaces[k].sourceDomain]->D2;
            switch (interfaces[k].sourceCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membranes[interfaces[k].sourceDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_south;
                    Zsource(k) = interfaces[k].lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interfaces[k].sourceDomain]->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interfaces[k].sourceDomain]->h_1s2_east;
                    Zsource(k) = (interfaces[k].lambdaSource*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interfaces[k].sourceDomain]->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membranes[interfaces[k].sourceDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_north;
                    Zsource(k) = interfaces[k].lambdaSource*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interfaces[k].sourceDomain]->h_1s2_west;
                    Zsource(k) = (interfaces[k].lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    break;
                }
            }
            Z  = membranes[interfaces[k].targetDomain]->z;
            D1 = membranes[interfaces[k].targetDomain]->D1;
            D2 = membranes[interfaces[k].targetDomain]->D2;
            switch (interfaces[k].targetCurve)
            {
                case 0: // South
                {
                    arma::vec dZd1 = D1*Z.col(0);
                    arma::vec dZd2 = Z*D2.row(0).t();
                    arma::vec h_2s1 = membranes[interfaces[k].targetDomain]->h_2s1_south;
                    arma::vec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_south;
                    Ztarget(k) = interfaces[k].lambdaSource*Z.col(0) - (h_2s1%dZd1 + h_2s2%dZd2);
                    if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 1: // East
                {
                    size_t nx = membranes[interfaces[k].targetDomain]->nx;
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_east;
                    arma::rowvec h_1s2 = membranes[interfaces[k].targetDomain]->h_1s2_east;
                    Ztarget(k) = (interfaces[k].lambdaSource*Z.row(nx-1) + h_1s1%dZd1 + h_1s2%dZd2).t();
                    if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 2: // North
                {
                    size_t ny = membranes[interfaces[k].targetDomain]->ny;
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membranes[interfaces[k].targetDomain]->h_2s1_north;
                    arma::vec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_north;
                    Ztarget(k) = interfaces[k].lambdaSource*Z.col(ny-1) + h_2s1%dZd1 + h_2s2%dZd2;
                    if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec dZd1 = D1.row(0)*Z;
                    arma::rowvec dZd2 = Z.row(0)*D2.t();
                    arma::rowvec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_west;
                    arma::rowvec h_1s2 = membranes[interfaces[k].targetDomain]->h_1s2_west;
                    Ztarget(k) = (interfaces[k].lambdaSource*Z.row(0) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
                        Ztarget(k) = reverse(Ztarget(k));
                    break;
                }
            }
        }
        converged = true;
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            printf("Interface %lu: ", k+1);
            double res = fabs(1 - arma::norm(Zsource(k).subvec(1, Zsource(k).size()-2))/arma::norm(Ztarget(k).subvec(1, Ztarget(k).size()-2)));
            printf("Residual %4.2e\n", res);
            if (res > residualTarget)
                converged = false;
        }
        count++;
        if (count > iterations)
            converged = true;
    } while (converged == false);
}