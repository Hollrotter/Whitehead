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
        std::cout << "Substep " << substep << '/' << substeps << '\n';
        for (size_t k = 0; k < membranes.size(); k++)
        {
            membranes[k]->p  =  p0(k)*substep/substeps;
            membranes[k]->p1 = p10(k)*substep/substeps;
            membranes[k]->p2 = p20(k)*substep/substeps;
        }
        bool converged = false;
        size_t count = 1;
        do
        {
            std::cout << "Iteration " << count << '/' << iterations << '\n';
            for (Interface& interface:interfaces)
            {
                Direction targetDirection = static_cast<Direction>(interface.targetCurve);
                Membrane *membraneSource = membranes[interface.sourceDomain];
                Membrane *membraneTarget = membranes[interface.targetDomain];
                arma::mat D1    = membraneSource->D1;
                arma::mat D2    = membraneSource->D2;
                arma::mat Z     = membraneSource->z;
                arma::mat V1    = membraneSource->v1;
                arma::mat V2    = membraneSource->v2;
                arma::mat v1__1 = membraneSource->v1__1;
                arma::mat v1__2 = membraneSource->v1__2;
                arma::mat v2__1 = membraneSource->v2__1;
                arma::mat v2__2 = membraneSource->v2__2;
                switch (interface.sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membraneSource->h_1s1_south;
                        arma::vec h_2s1 = membraneSource->h_2s1_south;
                        arma::vec h_2s2 = membraneSource->h_2s2_south;
                        arma::vec sourceZ  = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(0)
                                           + interface.c1Source%v1__1.col(0) + interface.c2Source%v1__2.col(0);
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                           + interface.c11Source%v1__1.col(0) + interface.c22Source%v2__2.col(0) + interface.c12Source%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                break;
                            case 1: // East
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 2: // North
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                break;
                            case 3: // West
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, sourceV2);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
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
                        arma::rowvec h_2s2 = membraneSource->h_2s2_east;
                        arma::vec sourceZ  = (interface.lambdaSource*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                           - interface.c11Source%v1__1.row(nx-1).t() - interface.c22Source%v2__2.row(nx-1).t() - interface.c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                           - interface.c1Source%v2__1.row(nx-1).t() - interface.c2Source%v2__2.row(nx-1).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 1: // East
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 2: // North
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, sourceV1);
                                break;
                            case 3: // West
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membraneSource->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membraneSource->h_1s1_north;
                        arma::vec h_2s1 = membraneSource->h_2s1_north;
                        arma::vec h_2s2 = membraneSource->h_2s2_north;
                        arma::vec sourceZ  = interface.lambdaSource*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                           - interface.c1Source%v1__1.col(ny-1) - interface.c2Source%v1__2.col(ny-1);
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                           - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV1);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                break;
                            case 1: // East
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, sourceV2);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                break;
                            case 2: // North
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV1)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 3: // West
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
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
                        arma::rowvec h_2s2 = membraneSource->h_2s2_west;
                        arma::vec sourceZ  = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                           + interface.c11Target%v1__1.row(0).t() + interface.c22Source%v2__2.row(0).t() + interface.c12Source%(v1__2.row(0) + v2__1.row(0)).t();
                        arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                           + interface.c1Source%v2__1.row(0).t() + interface.c2Source%v2__2.row(0).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1, sourceV2);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, sourceV1);
                                break;
                            case 1: // East
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, sourceZ);
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV1);
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, sourceV2);
                                break;
                            case 2: // North
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1, arma::vec(reverse(sourceV2)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                break;
                            case 3: // West
                                membraneTarget->boundary(Field::z,  targetDirection, BC::Robin, interface.lambdaSource,-1, arma::vec(reverse(sourceZ)));
                                membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV1)));
                                membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1, arma::vec(reverse(sourceV2)));
                                break;
                        }
                        break;
                    }
                }
                Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
                D1    = membraneTarget->D1;
                D2    = membraneTarget->D2;
                Z     = membraneTarget->z;
                V1    = membraneTarget->v1;
                V2    = membraneTarget->v2;
                v1__1 = membraneTarget->v1__1;
                v1__2 = membraneTarget->v1__2;
                v2__1 = membraneTarget->v2__1;
                v2__2 = membraneTarget->v2__2;
                switch (interface.targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membraneTarget->h_1s1_south;
                        arma::vec h_2s1 = membraneTarget->h_2s1_south;
                        arma::vec h_2s2 = membraneTarget->h_2s2_south;
                        arma::vec targetZ  = interface.lambdaTarget*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(0)
                                           + interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0);
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                           + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                break;
                            case 1: // East
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 2: // North
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                break;
                            case 3: // West
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, targetV2);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
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
                        arma::rowvec h_2s2 = membraneTarget->h_2s2_east;
                        arma::vec targetZ  = (interface.lambdaTarget*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                           - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(nx-1)).t()
                                           - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 1: // East
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 2: // North
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, targetV1);
                                break;
                            case 3: // West
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membraneTarget->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membraneTarget->h_1s1_north;
                        arma::vec h_2s1 = membraneTarget->h_2s1_north;
                        arma::vec h_2s2 = membraneTarget->h_2s2_north;
                        arma::vec targetZ  = interface.lambdaTarget*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(ny-1)
                                           - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                           - interface.c11Target%v1__1.col(ny-1) - interface.c22Target%v2__2.col(ny-1) - interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV1);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                break;
                            case 1: // East
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, targetV2);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                break;
                            case 2: // North
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV1)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 3: // West
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
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
                        arma::rowvec h_2s2 = membraneTarget->h_2s2_west;
                        arma::vec targetZ  = (interface.lambdaTarget*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                           + interface.c11Target%v1__1.row(0).t() + interface.c22Target%v2__2.row(0).t() + interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                        arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(0)).t()
                                           + interface.c1Target%v2__1.row(0).t() + interface.c2Target%v2__2.row(0).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1, targetV2);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, targetV1);
                                break;
                            case 1: // East
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetZ);
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV1);
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, targetV2);
                                break;
                            case 2: // North
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1, arma::vec(reverse(targetV2)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                break;
                            case 3: // West
                                membraneSource->boundary(Field::z,  sourceDirection, BC::Robin, interface.lambdaTarget,-1, arma::vec(reverse(targetZ)));
                                membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV1)));
                                membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1, arma::vec(reverse(targetV2)));
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
                Interface interface = interfaces[k];
                const Membrane *membraneSource = membranes[interface.sourceDomain];
                arma::mat D1    = membraneSource->D1;
                arma::mat D2    = membraneSource->D2;
                arma::mat Z     = membraneSource->z;
                arma::mat V1    = membraneSource->v1;
                arma::mat V2    = membraneSource->v2;
                arma::mat v1__1 = membraneSource->v1__1;
                arma::mat v1__2 = membraneSource->v1__2;
                arma::mat v2__1 = membraneSource->v2__1;
                arma::mat v2__2 = membraneSource->v2__2;
                switch (interface.sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membraneSource->h_1s1_south;
                        arma::vec h_2s1 = membraneSource->h_2s1_south;
                        arma::vec h_2s2 = membraneSource->h_2s2_south;
                        Zsource(k)  = interface.lambdaSource*Z.col(0) + h_2s1%dZd1 + h_2s2%dZd2;
                        Vtsource(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                    + interface.c1Source%v1__1.col(0) + interface.c2Source%v1__2.col(0);
                        Vnsource(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    + interface.c11Source%v1__1.col(0) + interface.c22Source%v2__2.col(0) + interface.c12Source%(v1__2.col(0) + v2__1.col(0));
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membraneSource->nx;
                        arma::rowvec dZd1 = D1.row(nx-1)*Z;
                        arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                        arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                        arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                        arma::rowvec h_2s2 = membraneSource->h_2s2_east;
                        Zsource(k)  = (interface.lambdaSource*Z.row(nx-1) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        Vnsource(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    - interface.c11Source%v1__1.row(nx-1).t() - interface.c22Source%v2__2.row(nx-1).t() - interface.c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                        Vtsource(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    - interface.c1Source%v2__1.row(nx-1).t() - interface.c2Source%v2__2.row(nx-1).t();
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membraneSource->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membraneSource->h_1s1_north;
                        arma::vec h_2s1 = membraneSource->h_2s1_north;
                        arma::vec h_2s2 = membraneSource->h_2s2_north;
                        Zsource(k)  = interface.lambdaSource*Z.col(ny-1) - h_2s1%dZd1 - h_2s2%dZd2;
                        Vtsource(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                    - interface.c1Source%v1__1.col(ny-1) - interface.c2Source%v1__2.col(ny-1);
                        Vnsource(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec dZd1 = D1.row(0)*Z;
                        arma::rowvec dZd2 = Z.row(0)*D2.t();
                        arma::rowvec h_1s1 = membraneSource->h_1s1_west;
                        arma::rowvec h_1s2 = membraneSource->h_1s2_west;
                        arma::rowvec h_2s2 = membraneSource->h_2s2_west;
                        Zsource(k)  = (interface.lambdaSource*Z.row(0) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        Vnsource(k) = (interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                    + interface.c11Source%v1__1.row(0).t() + interface.c22Source%v2__2.row(0).t() + interface.c12Source%(v1__2.row(0) + v2__1.row(0)).t()).t();
                        Vtsource(k) = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                    + interface.c1Source%v2__1.row(0).t() + interface.c2Source%v2__2.row(0).t();
                        break;
                    }
                }
                const Membrane *membraneTarget = membranes[interface.targetDomain];
                D1    = membraneTarget->D1;
                D2    = membraneTarget->D2;
                Z     = membraneTarget->z;
                V1    = membraneTarget->v1;
                V2    = membraneTarget->v2;
                v1__1 = membraneTarget->v1__1;
                v1__2 = membraneTarget->v1__2;
                v2__1 = membraneTarget->v2__1;
                v2__2 = membraneTarget->v2__2;
                switch (interface.targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec dZd1 = D1*Z.col(0);
                        arma::vec dZd2 = Z*D2.row(0).t();
                        arma::vec h_1s1 = membraneTarget->h_1s1_south;
                        arma::vec h_2s1 = membraneTarget->h_2s1_south;
                        arma::vec h_2s2 = membraneTarget->h_2s2_south;
                        Ztarget(k)  = interface.lambdaSource*Z.col(0) - h_2s1%dZd1 - h_2s2%dZd2;
                        if (interface.sourceCurve == 0) // South
                        {
                            Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(0)
                                        + interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0);
                            Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else if (interface.sourceCurve == 1) // East
                        {
                            Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(0)
                                        + interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0);
                            Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        - interface.c11Target%v1__1.col(0) - interface.c22Target%v2__2.col(0) - interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else if (interface.sourceCurve == 2) // North
                        {
                            Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                        - interface.c1Target%v1__1.col(0) - interface.c2Target%v1__2.col(0);
                            Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        - interface.c11Target%v1__1.col(0) - interface.c22Target%v2__2.col(0) - interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        else // West
                        {
                            Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                        - interface.c1Target%v1__1.col(0) - interface.c2Target%v1__2.col(0);
                            Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                        + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                        }
                        if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
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
                        arma::rowvec h_2s2 = membraneTarget->h_2s2_east;
                        Ztarget(k)  = (interface.lambdaSource*Z.row(nx-1) + h_1s1%dZd1 + h_1s2%dZd2).t();
                        if (interface.sourceCurve == 0) // South
                        {
                            Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        + interface.c11Target%v1__1.row(nx-1).t() + interface.c22Target%v2__2.row(nx-1).t() + interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        }
                        else if (interface.sourceCurve == 1) // East
                        {
                            Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        }
                        else if (interface.sourceCurve == 2) // North
                        {
                            Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        + interface.c1Target%v2__1.row(nx-1).t() + interface.c2Target%v2__2.row(nx-1).t();
                        }
                        else // West
                        {
                            Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                        + interface.c11Target%v1__1.row(nx-1).t() + interface.c22Target%v2__2.row(nx-1).t() + interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                            Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                        + interface.c1Target%v2__1.row(nx-1).t() + interface.c2Target%v2__2.row(nx-1).t();
                        }
                        if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                        {
                            Ztarget(k)  = reverse(Ztarget(k));
                            Vttarget(k) = reverse(Vttarget(k));
                            Vntarget(k) = reverse(Vntarget(k));
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membraneTarget->ny;
                        arma::vec dZd1 = D1*Z.col(ny-1);
                        arma::vec dZd2 = Z*D2.row(ny-1).t();
                        arma::vec h_1s1 = membraneTarget->h_1s1_north;
                        arma::vec h_2s1 = membraneTarget->h_2s1_north;
                        arma::vec h_2s2 = membraneTarget->h_2s2_north;
                        Ztarget(k)  = interface.lambdaSource*Z.col(ny-1) + h_2s1%dZd1 + h_2s2%dZd2;
                        if (interface.sourceCurve == 0) // South
                        {
                            Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                        + interface.c1Target%v1__1.col(ny-1) + interface.c2Target%v1__2.col(ny-1);
                            Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        + interface.c11Source%v1__1.col(ny-1) + interface.c22Source%v2__2.col(ny-1) + interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else if (interface.sourceCurve == 1) // East
                        {
                            Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                        + interface.c1Target%v1__1.col(ny-1) + interface.c2Target%v1__2.col(ny-1);
                            Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else if (interface.sourceCurve == 2) // North
                        {
                            Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(ny-1)
                                        - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                            Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        else // West
                        {
                            Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(ny-1)
                                        - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                            Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                        + interface.c11Source%v1__1.col(ny-1) + interface.c22Source%v2__2.col(ny-1) + interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        }
                        if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
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
                        arma::rowvec h_1s1 = membraneTarget->h_1s1_west;
                        arma::rowvec h_1s2 = membraneTarget->h_1s2_west;
                        arma::rowvec h_2s2 = membraneTarget->h_2s2_west;
                        Ztarget(k)  = (interface.lambdaSource*Z.row(0) - h_1s1%dZd1 - h_1s2%dZd2).t();
                        if (interface.sourceCurve == 0) // South
                        {
                            Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        + interface.c11Target%v1__1.row(0).t() + interface.c22Target%v2__2.row(0).t() + interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                        - interface.c1Target%v2__1.row(0).t() - interface.c2Target%v2__2.row(0).t();
                        }
                        else if (interface.sourceCurve == 1) // East
                        {
                            Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        - interface.c11Target%v1__1.row(0).t() - interface.c22Target%v2__2.row(0).t() - interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                        - interface.c1Target%v2__1.row(0).t() - interface.c2Target%v2__2.row(0).t();
                        }
                        else if (interface.sourceCurve == 2) // North
                        {
                            Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        - interface.c11Target%v1__1.row(0).t() - interface.c22Target%v2__2.row(0).t() - interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                        + interface.c1Target%v2__1.row(0).t() + interface.c2Target%v2__2.row(0).t();
                        }
                        else // West
                        {
                            Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                        + interface.c11Target%v1__1.row(0).t() + interface.c22Target%v2__2.row(0).t() + interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                            Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                        + interface.c1Target%v2__1.row(0).t() + interface.c2Target%v2__2.row(0).t();
                        }
                        if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
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
                std::cout << "Interface " << k+1 << '\n';
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