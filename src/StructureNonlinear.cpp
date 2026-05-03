#include "Structure.hpp"

void Structure::nonlinear()
{
    analysis = Analysis::nonlinear;
    arma::vec iter_old(membranes.size());
    arma::field<arma::mat>  p0(membranes.size());
    arma::field<arma::mat> p10(membranes.size());
    arma::field<arma::mat> p20(membranes.size());
    for (size_t k = 0; k < membranes.size(); k++)
    {
        iter_old(k) = membranes[k]->iter;
        membranes[k]->iter = 1;
        p0(k)  = membranes[k]->p;
        p10(k) = membranes[k]->p1;
        p20(k) = membranes[k]->p2;
    }
    arma::field<arma::vec> Vntarget(interfaces.size()), Vnsource(interfaces.size());
    arma::field<arma::vec> Vttarget(interfaces.size()), Vtsource(interfaces.size());
    arma::field<arma::vec>  Ztarget(interfaces.size()),  Zsource(interfaces.size());

    for (size_t substep = 1; substep <= substeps; substep++)
    {
        printf("Substep %lu/%lu\n", substep, substeps);
        for (size_t k = 0; k < membranes.size(); k++)
        {
            membranes[k]->p  =  p0(k)*substep/substeps;
            membranes[k]->p1 = p10(k)*substep/substeps;
            membranes[k]->p2 = p20(k)*substep/substeps;
        }
        bool converged = false;
        int count = 1;
        do
        {
            std::cout << "Iteration " << count << '/' << iterations << '\n';
            for (Interface& interface:interfaces)
            {
                Direction targetDirection = static_cast<Direction>(interface.targetCurve);
                Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
                arma::mat D1    = membranes[interface.sourceDomain]->D1;
                arma::mat D2    = membranes[interface.sourceDomain]->D2;
                arma::mat Z     = membranes[interface.sourceDomain]->z;
                arma::mat V1    = membranes[interface.sourceDomain]->v1;
                arma::mat V2    = membranes[interface.sourceDomain]->v2;
                arma::mat v1__1 = membranes[interface.sourceDomain]->v1__1;
                arma::mat v1__2 = membranes[interface.sourceDomain]->v1__2;
                arma::mat v2__1 = membranes[interface.sourceDomain]->v2__1;
                arma::mat v2__2 = membranes[interface.sourceDomain]->v2__2;
                switch (interface.sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membranes[interface.sourceDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_south;
                        arma::vec sourceZ  = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(0)
                                           + interface.c1Source%v1__1.col(0) + interface.c2Source%v1__2.col(0);
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                           + interface.c11Source%v1__1.col(0) + interface.c22Source%v2__2.col(0) + interface.c12Source%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
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
                        arma::rowvec h_2s2 = membranes[interface.sourceDomain]->h_2s2_east;
                        arma::vec sourceZ  = (interface.lambdaSource*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                           - interface.c11Source%v1__1.row(nx-1).t() - interface.c22Source%v2__2.row(nx-1).t() - interface.c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                           - interface.c1Source%v2__1.row(nx-1).t() - interface.c2Source%v2__2.row(nx-1).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, sourceV1);
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interface.sourceDomain]->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membranes[interface.sourceDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_north;
                        arma::vec sourceZ  = interface.lambdaSource*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                           - interface.c1Source%v1__1.col(ny-1) - interface.c2Source%v1__2.col(ny-1);
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                           - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
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
                        arma::rowvec h_2s2 = membranes[interface.sourceDomain]->h_2s2_west;
                        arma::vec sourceZ  = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                           + interface.c11Target%v1__1.row(0).t() + interface.c22Source%v2__2.row(0).t() + interface.c12Source%(v1__2.row(0) + v2__1.row(0)).t();
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                           + interface.c1Source%v2__1.row(0).t() + interface.c2Source%v2__2.row(0).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, sourceV1);
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                break;
                        }
                        break;
                    }
                }
                D1    = membranes[interface.targetDomain]->D1;
                D2    = membranes[interface.targetDomain]->D2;
                Z     = membranes[interface.targetDomain]->z;
                V1    = membranes[interface.targetDomain]->v1;
                V2    = membranes[interface.targetDomain]->v2;
                v1__1 = membranes[interface.targetDomain]->v1__1;
                v1__2 = membranes[interface.targetDomain]->v1__2;
                v2__1 = membranes[interface.targetDomain]->v2__1;
                v2__2 = membranes[interface.targetDomain]->v2__2;
                switch (interface.targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membranes[interface.targetDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_south;
                        arma::vec targetZ  = interface.lambdaTarget*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(0)
                                           + interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0);
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                           + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
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
                        arma::rowvec h_2s2 = membranes[interface.targetDomain]->h_2s2_east;
                        arma::vec targetZ  = (interface.lambdaTarget*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                           - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(nx-1)).t()
                                           - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, targetV1);
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interface.targetDomain]->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membranes[interface.targetDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_north;
                        arma::vec targetZ  = interface.lambdaTarget*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(ny-1)
                                           - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                           - interface.c11Target%v1__1.col(ny-1) - interface.c22Target%v2__2.col(ny-1) - interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
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
                        arma::rowvec h_2s2 = membranes[interface.targetDomain]->h_2s2_west;
                        arma::vec targetZ  = (interface.lambdaTarget*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                           + interface.c11Target%v1__1.row(0).t() + interface.c22Target%v2__2.row(0).t() + interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(0)).t()
                                           + interface.c1Target%v2__1.row(0).t() + interface.c2Target%v2__2.row(0).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, targetV1);
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                break;
                        }
                        break;
                    }
                }
            }
            #pragma omp parallel for
            for (auto& membrane:membranes)
                membrane->nonlinear();
            
            for (size_t k = 0; k < interfaces.size(); k++)
            {
                arma::mat D1    = membranes[interfaces[k].sourceDomain]->D1;
                arma::mat D2    = membranes[interfaces[k].sourceDomain]->D2;
                arma::mat Z     = membranes[interfaces[k].sourceDomain]->z;
                arma::mat V1    = membranes[interfaces[k].sourceDomain]->v1;
                arma::mat V2    = membranes[interfaces[k].sourceDomain]->v2;
                arma::mat v1__1 = membranes[interfaces[k].sourceDomain]->v1__1;
                arma::mat v1__2 = membranes[interfaces[k].sourceDomain]->v1__2;
                arma::mat v2__1 = membranes[interfaces[k].sourceDomain]->v2__1;
                arma::mat v2__2 = membranes[interfaces[k].sourceDomain]->v2__2;
                switch (interfaces[k].sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interfaces[k].sourceDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_south;
                        Zsource(k)  = interfaces[k].lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        Vtsource(k) = interfaces[k].lambdaSource*h_1s1%V1.col(0)
                                    + interfaces[k].c1Source%v1__1.col(0) + interfaces[k].c2Source%v1__2.col(0);
                        Vnsource(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    + interfaces[k].c11Source%v1__1.col(0) + interfaces[k].c22Source%v2__2.col(0) + interfaces[k].c12Source%(v1__2.col(0) + v2__1.col(0));
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interfaces[k].sourceDomain]->nx;
                        arma::rowvec dZd1 = D1.row(nx-1)*Z;
                        arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                        arma::rowvec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_east;
                        arma::rowvec h_1s2 = membranes[interfaces[k].sourceDomain]->h_1s2_east;
                        arma::rowvec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_east;
                        Zsource(k)  = (interfaces[k].lambdaSource*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        Vnsource(k) = interfaces[k].lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    - interfaces[k].c11Source%v1__1.row(nx-1).t() - interfaces[k].c22Source%v2__2.row(nx-1).t() - interfaces[k].c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        Vtsource(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    - interfaces[k].c1Source%v2__1.row(nx-1).t() - interfaces[k].c2Source%v2__2.row(nx-1).t();
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interfaces[k].sourceDomain]->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interfaces[k].sourceDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_north;
                        Zsource(k)  = interfaces[k].lambdaSource*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        Vtsource(k) = interfaces[k].lambdaSource*h_1s1%V1.col(ny-1)
                                    - interfaces[k].c1Source%v1__1.col(ny-1) - interfaces[k].c2Source%v1__2.col(ny-1);
                        Vnsource(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    - interfaces[k].c11Source%v1__1.col(ny-1) - interfaces[k].c22Source%v2__2.col(ny-1) - interfaces[k].c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec dZd1 = D1.row(0)*Z;
                        arma::rowvec dZd2 = Z.row(0)*D2.t();
                        arma::rowvec h_1s1 = membranes[interfaces[k].sourceDomain]->h_1s1_west;
                        arma::rowvec h_1s2 = membranes[interfaces[k].sourceDomain]->h_1s2_west;
                        arma::rowvec h_2s2 = membranes[interfaces[k].sourceDomain]->h_2s2_west;
                        Zsource(k)  = (interfaces[k].lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        Vnsource(k) = (interfaces[k].lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                    + interfaces[k].c11Source%v1__1.row(0).t() + interfaces[k].c22Source%v2__2.row(0).t() + interfaces[k].c12Source%(v1__2.row(0) + v2__1.row(0)).t()).t();
                        Vtsource(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(0)).t()
                                    + interfaces[k].c1Source%v2__1.row(0).t() + interfaces[k].c2Source%v2__2.row(0).t();
                        break;
                    }
                }
                D1    = membranes[interfaces[k].targetDomain]->D1;
                D2    = membranes[interfaces[k].targetDomain]->D2;
                Z     = membranes[interfaces[k].targetDomain]->z;
                V1    = membranes[interfaces[k].targetDomain]->v1;
                V2    = membranes[interfaces[k].targetDomain]->v2;
                v1__1 = membranes[interfaces[k].targetDomain]->v1__1;
                v1__2 = membranes[interfaces[k].targetDomain]->v1__2;
                v2__1 = membranes[interfaces[k].targetDomain]->v2__1;
                v2__2 = membranes[interfaces[k].targetDomain]->v2__2;
                switch (interfaces[k].targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interfaces[k].targetDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_south;
                        Ztarget(k)  = interfaces[k].lambdaSource*Z.col(0) - h_2s1%dZd1 - h_2s2%dZd2;
                        if (interfaces[k].sourceCurve == 0) // South
                        {
                            Vttarget(k) =-interfaces[k].lambdaSource*h_1s1%V1.col(0)
                                        + interfaces[k].c1Target%v1__1.col(0) + interfaces[k].c2Target%v1__2.col(0);
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        + interfaces[k].c11Target%v1__1.col(0) + interfaces[k].c22Target%v2__2.col(0) + interfaces[k].c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else if (interfaces[k].sourceCurve == 1) // East
                        {
                            Vttarget(k) =-interfaces[k].lambdaSource*h_1s1%V1.col(0)
                                        + interfaces[k].c1Target%v1__1.col(0) + interfaces[k].c2Target%v1__2.col(0);
                            Vntarget(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        - interfaces[k].c11Target%v1__1.col(0) - interfaces[k].c22Target%v2__2.col(0) - interfaces[k].c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else if (interfaces[k].sourceCurve == 2) // North
                        {
                            Vttarget(k) = interfaces[k].lambdaSource*h_1s1%V1.col(0)
                                        - interfaces[k].c1Target%v1__1.col(0) - interfaces[k].c2Target%v1__2.col(0);
                            Vntarget(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        - interfaces[k].c11Target%v1__1.col(0) - interfaces[k].c22Target%v2__2.col(0) - interfaces[k].c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else // West
                        {
                            Vttarget(k) = interfaces[k].lambdaSource*h_1s1%V1.col(0)
                                        - interfaces[k].c1Target%v1__1.col(0) - interfaces[k].c2Target%v1__2.col(0);
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        + interfaces[k].c11Target%v1__1.col(0) + interfaces[k].c22Target%v2__2.col(0) + interfaces[k].c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
                        }
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interfaces[k].targetDomain]->nx;
                        arma::rowvec dZd1 = D1.row(nx-1)*Z;
                        arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                        arma::rowvec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_east;
                        arma::rowvec h_1s2 = membranes[interfaces[k].targetDomain]->h_1s2_east;
                        arma::rowvec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_east;
                        Ztarget(k)  = (interfaces[k].lambdaSource*Z.row(nx-1) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        if (interfaces[k].sourceCurve == 0) // South
                        {
                            Vntarget(k) = interfaces[k].lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        + interfaces[k].c11Target%v1__1.row(nx-1).t() + interfaces[k].c22Target%v2__2.row(nx-1).t() + interfaces[k].c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) =-interfaces[k].lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        - interfaces[k].c1Target%v2__1.row(nx-1).t() - interfaces[k].c2Target%v2__2.row(nx-1).t();
                        }
                        else if (interfaces[k].sourceCurve == 1) // East
                        {
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        - interfaces[k].c11Target%v1__1.row(nx-1).t() - interfaces[k].c22Target%v2__2.row(nx-1).t() - interfaces[k].c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) =-interfaces[k].lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        - interfaces[k].c1Target%v2__1.row(nx-1).t() - interfaces[k].c2Target%v2__2.row(nx-1).t();
                        }
                        else if (interfaces[k].sourceCurve == 2) // North
                        {
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        - interfaces[k].c11Target%v1__1.row(nx-1).t() - interfaces[k].c22Target%v2__2.row(nx-1).t() - interfaces[k].c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        + interfaces[k].c1Target%v2__1.row(nx-1).t() + interfaces[k].c2Target%v2__2.row(nx-1).t();
                        }
                        else // West
                        {
                            Vntarget(k) = interfaces[k].lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        + interfaces[k].c11Target%v1__1.row(nx-1).t() + interfaces[k].c22Target%v2__2.row(nx-1).t() + interfaces[k].c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        + interfaces[k].c1Target%v2__1.row(nx-1).t() + interfaces[k].c2Target%v2__2.row(nx-1).t();
                        }
                        if (interfaces[k].sourceCurve == 0 || interfaces[k].sourceCurve == 1)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interfaces[k].targetDomain]->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interfaces[k].targetDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_north;
                        Ztarget(k)  = interfaces[k].lambdaSource*Z.col(ny-1) + h_2s1%dZd1 + h_2s2%dZd2;
                        if (interfaces[k].sourceCurve == 0) // South
                        {
                            Vttarget(k) = interfaces[k].lambdaSource*h_1s1%V1.col(ny-1)
                                        + interfaces[k].c1Target%v1__1.col(ny-1) + interfaces[k].c2Target%v1__2.col(ny-1);
                            Vntarget(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        + interfaces[k].c11Source%v1__1.col(ny-1) + interfaces[k].c22Source%v2__2.col(ny-1) + interfaces[k].c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else if (interfaces[k].sourceCurve == 1) // East
                        {
                            Vttarget(k) = interfaces[k].lambdaSource*h_1s1%V1.col(ny-1)
                                        + interfaces[k].c1Target%v1__1.col(ny-1) + interfaces[k].c2Target%v1__2.col(ny-1);
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        - interfaces[k].c11Source%v1__1.col(ny-1) - interfaces[k].c22Source%v2__2.col(ny-1) - interfaces[k].c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else if (interfaces[k].sourceCurve == 2) // North
                        {
                            Vttarget(k) =-interfaces[k].lambdaSource*h_1s1%V1.col(ny-1)
                                        - interfaces[k].c1Target%v1__1.col(ny-1) - interfaces[k].c2Target%v1__2.col(ny-1);
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        - interfaces[k].c11Source%v1__1.col(ny-1) - interfaces[k].c22Source%v2__2.col(ny-1) - interfaces[k].c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else // West
                        {
                            Vttarget(k) =-interfaces[k].lambdaSource*h_1s1%V1.col(ny-1)
                                        - interfaces[k].c1Target%v1__1.col(ny-1) - interfaces[k].c2Target%v1__2.col(ny-1);
                            Vntarget(k) = interfaces[k].lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        + interfaces[k].c11Source%v1__1.col(ny-1) + interfaces[k].c22Source%v2__2.col(ny-1) + interfaces[k].c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
                        }
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec dZd1 = D1.row(0)*Z;
                        arma::rowvec dZd2 = Z.row(0)*D2.t();
                        arma::rowvec h_1s1 = membranes[interfaces[k].targetDomain]->h_1s1_west;
                        arma::rowvec h_1s2 = membranes[interfaces[k].targetDomain]->h_1s2_west;
                        arma::rowvec h_2s2 = membranes[interfaces[k].targetDomain]->h_2s2_west;
                        Ztarget(k)  = (interfaces[k].lambdaSource*Z.row(0) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        if (interfaces[k].sourceCurve == 0) // South
                        {
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        + interfaces[k].c11Target%v1__1.row(0).t() + interfaces[k].c22Target%v2__2.row(0).t() + interfaces[k].c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(0)).t()
                                        - interfaces[k].c1Target%v2__1.row(0).t() - interfaces[k].c2Target%v2__2.row(0).t();
                        }
                        else if (interfaces[k].sourceCurve == 1) // East
                        {
                            Vntarget(k) = interfaces[k].lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        - interfaces[k].c11Target%v1__1.row(0).t() - interfaces[k].c22Target%v2__2.row(0).t() - interfaces[k].c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) = interfaces[k].lambdaSource*(h_2s2%V2.row(0)).t()
                                        - interfaces[k].c1Target%v2__1.row(0).t() - interfaces[k].c2Target%v2__2.row(0).t();
                        }
                        else if (interfaces[k].sourceCurve == 2) // North
                        {
                            Vntarget(k) = interfaces[k].lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        - interfaces[k].c11Target%v1__1.row(0).t() - interfaces[k].c22Target%v2__2.row(0).t() - interfaces[k].c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) =-interfaces[k].lambdaSource*(h_2s2%V2.row(0)).t()
                                        + interfaces[k].c1Target%v2__1.row(0).t() + interfaces[k].c2Target%v2__2.row(0).t();
                        }
                        else // West
                        {
                            Vntarget(k) =-interfaces[k].lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        + interfaces[k].c11Target%v1__1.row(0).t() + interfaces[k].c22Target%v2__2.row(0).t() + interfaces[k].c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) =-interfaces[k].lambdaSource*(h_2s2%V2.row(0)).t()
                                        + interfaces[k].c1Target%v2__1.row(0).t() + interfaces[k].c2Target%v2__2.row(0).t();
                        }
                        if (interfaces[k].sourceCurve == 2 || interfaces[k].sourceCurve == 3)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
                        }
                        break;
                    }
                }
            }
            converged = true;
            for (size_t k = 0; k < interfaces.size(); k++)
            {
                printf("Interface %lu\n", k+1);
                double res1 = fabs(1 - arma::norm(Vnsource(k).subvec(1, Vnsource(k).size()-2))/arma::norm(Vntarget(k).subvec(1, Vntarget(k).size()-2)));
                double res2 = fabs(1 - arma::norm(Vtsource(k).subvec(1, Vtsource(k).size()-2))/arma::norm(Vttarget(k).subvec(1, Vttarget(k).size()-2)));
                double res3 = fabs(1 - arma::norm( Zsource(k).subvec(1,  Zsource(k).size()-2))/arma::norm( Ztarget(k).subvec(1,  Ztarget(k).size()-2)));
                printf("Residual normal  %4.2e\n", res1);
                printf("Residual tangent %4.2e\n", res2);
                printf("Residual z       %4.2e\n", res3);
                if (res1 > residualTarget || res2 > residualTarget || res3 > residualTarget)
                    converged = false;
            }
            count++;
            if (count > iterations)
                converged = true;
        } while (converged == false);
    }
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->iter = iter_old(k);
}