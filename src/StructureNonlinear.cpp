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
    arma::field<arma::vec> Ztarget(interfaces.size()),  Zsource(interfaces.size());
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
        int count = 0;
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
                arma::mat V1   = membranes[interface.sourceDomain]->v1;
                arma::mat V2   = membranes[interface.sourceDomain]->v2;
                arma::mat v1__1 = membranes[interface.sourceDomain]->v1__1;
                arma::mat v1__2 = membranes[interface.sourceDomain]->v1__2;
                arma::mat v2__1 = membranes[interface.sourceDomain]->v2__1;
                arma::mat v2__2 = membranes[interface.sourceDomain]->v2__2;
                arma::mat dZd1  = D1*Z;
                arma::mat dZd2  = Z*D2.t();
                double gamma;
                arma::vec sourceZ;
                switch (interface.sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec h_1s1 = membranes[interface.sourceDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_south;
                        gamma = gamma0*mean(h_2s2);
                        sourceZ = gamma*Z.col(0) + h_2s1%dZd1.col(0) + h_2s2%dZd2.col(0);
                        arma::vec c1 = h_1s1%h_2s1;
                        arma::vec c2 = h_1s1%h_2s2;
                        arma::vec sourceV1 = gamma*h_1s1%V1.col(0) + c1%v1__1.col(0) + c2%v1__2.col(0);
                        arma::vec c11 = pow(h_2s1, 2);
                        arma::vec c12 = h_2s1%h_2s2;
                        arma::vec c22 = pow(h_2s2, 2);
                        arma::vec sourceV2 = gamma*(h_2s1%V1.col(0) + h_2s2%V2.col(0)) + c11%v1__1.col(0) + c22%v2__2.col(0) + c12%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV2)));
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma, 1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma, 1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma, 1, sourceV2);
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma, 1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma,-1, sourceV1);
                                break;
                        }
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interface.sourceDomain]->nx;
                        arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_east;
                        arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_east;
                        arma::rowvec h_2s2 = membranes[interface.sourceDomain]->h_2s2_east;
                        gamma = gamma0*mean(h_1s1);
                        sourceZ = (gamma*Z.row(nx-1) - h_1s1%dZd1.row(nx-1) - h_1s2%dZd2.row(nx-1)).t();
                        arma::rowvec c11 = pow(h_1s1, 2);
                        arma::rowvec c12 = h_1s1%h_1s2;
                        arma::rowvec c22 = pow(h_1s2, 2);
                        arma::vec sourceV1 = (gamma*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)) - c11%v1__1.row(nx-1) - c22%v2__2.row(nx-1) - c12%(v1__2.row(nx-1) + v2__1.row(nx-1))).t();
                        arma::rowvec c1 = h_2s2%h_1s1;
                        arma::rowvec c2 = h_2s2%h_1s2;
                        arma::vec sourceV2 = (gamma*h_2s2%V2.row(nx-1) - c1%v2__1.row(nx-1) - c2%v2__2.row(nx-1)).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma,-1, arma::vec(reverse(sourceV1)));
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma, 1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma,-1, sourceV1);
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma,-1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma,-1, sourceV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interface.sourceDomain]->ny;
                        arma::vec h_1s1 = membranes[interface.sourceDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interface.sourceDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interface.sourceDomain]->h_2s2_north;
                        gamma = gamma0*mean(h_2s2);
                        sourceZ = gamma*Z.col(ny-1) - h_2s1%dZd1.col(ny-1) - h_2s2%dZd2.col(ny-1);
                        arma::vec c1 = h_1s1%h_2s1;
                        arma::vec c2 = h_1s1%h_2s2;
                        arma::vec sourceV1 = gamma*h_1s1%V1.col(ny-1) - c1%v1__1.col(ny-1) - c2%v1__2.col(ny-1);
                        arma::vec c11 = pow(h_2s1, 2);
                        arma::vec c12 = h_2s1%h_2s2;
                        arma::vec c22 = pow(h_2s2, 2);
                        arma::vec sourceV2 = gamma*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1)) - c11%v1__1.col(ny-1) - c22%v2__2.col(ny-1) - c12%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma,-1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma,-1, sourceV2);
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma,-1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma, 1, sourceV1);
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV2)));
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma,-1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV1)));
                                break;
                        }
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec h_1s1 = membranes[interface.sourceDomain]->h_1s1_west;
                        arma::rowvec h_1s2 = membranes[interface.sourceDomain]->h_1s2_west;
                        arma::rowvec h_2s2 = membranes[interface.sourceDomain]->h_2s2_west;
                        gamma = gamma0*mean(h_1s1);
                        sourceZ = (gamma*Z.row(0) + h_1s1%dZd1.row(0) + h_1s2%dZd2.row(0)).t();
                        arma::rowvec c11 = pow(h_1s1, 2);
                        arma::rowvec c12 = h_1s1%h_1s2;
                        arma::rowvec c22 = pow(h_1s2, 2);
                        arma::vec sourceV1 = (gamma*(h_1s1%V1.row(0) + h_1s2%V2.row(0)) + c11%v1__1.row(0) + c22%v2__2.row(0) + c12%(v1__2.row(0) + v2__1.row(0))).t();
                        arma::rowvec c1 = h_2s2%h_1s1;
                        arma::rowvec c2 = h_2s2%h_1s2;
                        arma::vec sourceV2 = (gamma*h_2s2%V2.row(0) + c1%v2__1.row(0) + c2%v2__2.row(0)).t();
                        switch (interface.targetCurve)
                        {
                            case 0: // South
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma,-1, sourceV2);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma, 1, sourceV1);
                                break;
                            case 1: // East
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin, gamma, 1, sourceV1);
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma, 1, sourceV2);
                                break;
                            case 2: // North
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma,-1, arma::vec(reverse(sourceV2)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin, gamma, 1, arma::vec(reverse(sourceV1)));
                                break;
                            case 3: // West
                                membranes[interface.targetDomain]->boundary(Field::v1, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV1)));
                                membranes[interface.targetDomain]->boundary(Field::v2, targetDirection, BC::Robin,-gamma, 1, arma::vec(reverse(sourceV2)));
                                break;
                        }
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
                V1   = membranes[interface.targetDomain]->v1;
                V2   = membranes[interface.targetDomain]->v2;
                v1__1 = membranes[interface.targetDomain]->v1__1;
                v1__2 = membranes[interface.targetDomain]->v1__2;
                v2__1 = membranes[interface.targetDomain]->v2__1;
                v2__2 = membranes[interface.targetDomain]->v2__2;
                dZd1  = D1*Z;
                dZd2  = Z*D2.t();
                arma::vec targetZ;
                switch (interface.targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec h_1s1 = membranes[interface.targetDomain]->h_1s1_south;
                        arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_south;
                        arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_south;
                        gamma = gamma0*mean(h_2s2);
                        targetZ = gamma*Z.col(0) + h_2s1%dZd1.col(0) + h_2s2%dZd2.col(0);
                        arma::vec c1 = h_1s1%h_2s1;
                        arma::vec c2 = h_1s1%h_2s2;
                        arma::vec targetV1 = gamma*h_1s1%V1.col(0) + c1%v1__1.col(0) + c2%v1__2.col(0);
                        arma::vec c11 = pow(h_2s1, 2);
                        arma::vec c12 = h_2s1%h_2s2;
                        arma::vec c22 = pow(h_2s2, 2);
                        arma::vec targetV2 = gamma*(h_2s1%V1.col(0) + h_2s2%V2.col(0)) + c11%v1__1.col(0) + c22%v2__2.col(0) + c12%(v1__2.col(0) + v2__1.col(0));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV2)));
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma, 1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma, 1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma, 1, targetV2);
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma, 1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma,-1, targetV1);
                                break;
                        }
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interface.targetDomain]->nx;
                        arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_east;
                        arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_east;
                        arma::rowvec h_2s2 = membranes[interface.targetDomain]->h_2s2_east;
                        gamma = gamma0*mean(h_1s1);
                        targetZ = (gamma*Z.row(nx-1) - h_1s1%dZd1.row(nx-1) - h_1s2%dZd2.row(nx-1)).t();
                        arma::rowvec c11 = pow(h_1s1, 2);
                        arma::rowvec c12 = h_1s1%h_1s2;
                        arma::rowvec c22 = pow(h_1s2, 2);
                        arma::vec targetV1 = (gamma*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)) - c11%v1__1.row(nx-1) - c22%v2__2.row(nx-1) - c12%(v1__2.row(nx-1) + v2__1.row(nx-1))).t();
                        arma::rowvec c1 = h_2s2%h_1s1;
                        arma::rowvec c2 = h_2s2%h_1s2;
                        arma::vec targetV2 = (gamma*h_2s2%V2.row(nx-1) - c1%v2__1.row(nx-1) - c2%v2__2.row(nx-1)).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma,-1, arma::vec(reverse(targetV1)));
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma, 1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma,-1, targetV1);
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma,-1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma,-1, targetV2);
                                break;
                        }
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interface.targetDomain]->ny;
                        arma::vec h_1s1 = membranes[interface.targetDomain]->h_1s1_north;
                        arma::vec h_2s1 = membranes[interface.targetDomain]->h_2s1_north;
                        arma::vec h_2s2 = membranes[interface.targetDomain]->h_2s2_north;
                        gamma = gamma0*mean(h_2s2);
                        targetZ = gamma*Z.col(ny-1) - h_2s1%dZd1.col(ny-1) - h_2s2%dZd2.col(ny-1);
                        arma::vec c1 = h_1s1%h_2s1;
                        arma::vec c2 = h_1s1%h_2s2;
                        arma::vec targetV1 = gamma*h_1s1%V1.col(ny-1) - c1%v1__1.col(ny-1) - c2%v1__2.col(ny-1);
                        arma::vec c11 = pow(h_2s1, 2);
                        arma::vec c12 = h_2s1%h_2s2;
                        arma::vec c22 = pow(h_2s2, 2);
                        arma::vec targetV2 = gamma*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1)) - c11%v1__1.col(ny-1) - c22%v2__2.col(ny-1) - c12%(v1__2.col(ny-1) + v2__1.col(ny-1));
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma,-1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma,-1, targetV2);
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma,-1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma, 1, targetV1);
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV2)));
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma,-1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV1)));
                                break;
                        }
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec h_1s1 = membranes[interface.targetDomain]->h_1s1_west;
                        arma::rowvec h_1s2 = membranes[interface.targetDomain]->h_1s2_west;
                        arma::rowvec h_2s2 = membranes[interface.targetDomain]->h_2s2_west;
                        gamma = gamma0*mean(h_1s1);
                        targetZ = (gamma*Z.row(0) + h_1s1%dZd1.row(0) + h_1s2%dZd2.row(0)).t();
                        arma::rowvec c11 = pow(h_1s1, 2);
                        arma::rowvec c12 = h_1s1%h_1s2;
                        arma::rowvec c22 = pow(h_1s2, 2);
                        arma::vec targetV1 = (gamma*(h_1s1%V1.row(0) + h_1s2%V2.row(0)) + c11%v1__1.row(0) + c22%v2__2.row(0) + c12%(v1__2.row(0) + v2__1.row(0))).t();
                        arma::rowvec c1 = h_2s2%h_1s1;
                        arma::rowvec c2 = h_2s2%h_1s2;
                        arma::vec targetV2 = (gamma*h_2s2%V2.row(0) + c1%v2__1.row(0) + c2%v2__2.row(0)).t();
                        switch (interface.sourceCurve)
                        {
                            case 0: // South
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma,-1, targetV2);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma, 1, targetV1);
                                break;
                            case 1: // East
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin, gamma, 1, targetV1);
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma, 1, targetV2);
                                break;
                            case 2: // North
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma,-1, arma::vec(reverse(targetV2)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin, gamma, 1, arma::vec(reverse(targetV1)));
                                break;
                            case 3: // West
                                membranes[interface.sourceDomain]->boundary(Field::v1, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV1)));
                                membranes[interface.sourceDomain]->boundary(Field::v2, sourceDirection, BC::Robin,-gamma, 1, arma::vec(reverse(targetV2)));
                                break;
                        }
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
                membrane->nonlinear();
            count++;
            for (size_t k = 0; k < interfaces.size(); k++)
            {
                arma::mat e11  = membranes[interfaces[k].targetDomain]->e11();
                arma::mat e12  = membranes[interfaces[k].targetDomain]->e12();
                arma::mat e22  = membranes[interfaces[k].targetDomain]->e22();
                arma::mat e_11 = membranes[interfaces[k].targetDomain]->e_11();
                arma::mat e_22 = membranes[interfaces[k].targetDomain]->e_22();
                arma::mat V1   = membranes[interfaces[k].targetDomain]->v1;
                arma::mat V2   = membranes[interfaces[k].targetDomain]->v2;
                switch (interfaces[k].targetCurve)
                {
                    case 0: // South
                    {
                        arma::vec h_1s1 = 1/sqrt(e_11.col(0));
                        arma::vec h_2s1 = e12.col(0)/sqrt(e22.col(0));
                        arma::vec h_2s2 = sqrt(e22.col(0));
                        Ztarget(k) = membranes[interfaces[k].targetDomain]->z.col(0);
                        Vttarget(k) = h_1s1%V1.col(0);
                        Vntarget(k) = h_2s1%V1.col(0) + h_2s2%V2.col(0);
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interfaces[k].targetDomain]->nx;
                        arma::rowvec h_1s1 = sqrt(e11.row(nx-1));
                        arma::rowvec h_1s2 = e12.row(nx-1)/sqrt(e11.row(nx-1));
                        arma::rowvec h_2s2 = 1/sqrt(e_22.row(nx-1));
                        Ztarget(k) = membranes[interfaces[k].targetDomain]->z.row(nx-1).t();
                        Vntarget(k) = (h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t();
                        Vttarget(k) = (h_2s2%V2.row(nx-1)).t();
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interfaces[k].targetDomain]->ny;
                        arma::vec h_1s1 = 1/sqrt(e_11.col(ny-1));
                        arma::vec h_2s1 = e12.col(ny-1)/sqrt(e22.col(ny-1));
                        arma::vec h_2s2 = sqrt(e22.col(ny-1));
                        Ztarget(k) = membranes[interfaces[k].targetDomain]->z.col(ny-1);
                        Vttarget(k) = h_1s1%V1.col(ny-1);
                        Vntarget(k) = h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1);
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec h_1s1 = sqrt(e11.row(0));
                        arma::rowvec h_1s2 = e12.row(0)/sqrt(e11.row(0));
                        arma::rowvec h_2s2 = 1/sqrt(e_22.row(0));
                        Ztarget(k) = membranes[interfaces[k].targetDomain]->z.row(0).t();
                        Vntarget(k) = (h_1s1%V1.row(0) + h_1s2%V2.row(0)).t();
                        Vttarget(k) = (h_2s2%V2.row(0)).t();
                        break;
                    }
                }
                e11  = membranes[interfaces[k].sourceDomain]->e11();
                e12  = membranes[interfaces[k].sourceDomain]->e12();
                e22  = membranes[interfaces[k].sourceDomain]->e22();
                e_11 = membranes[interfaces[k].sourceDomain]->e_11();
                e_22 = membranes[interfaces[k].sourceDomain]->e_22();
                V1   = membranes[interfaces[k].sourceDomain]->v1;
                V2   = membranes[interfaces[k].sourceDomain]->v2;
                switch (interfaces[k].sourceCurve)
                {
                    case 0: // South
                    {
                        arma::vec h_1s1 = 1/sqrt(e_11.col(0));
                        arma::vec h_2s1 = e12.col(0)/sqrt(e22.col(0));
                        arma::vec h_2s2 = sqrt(e22.col(0));
                        Zsource(k) = membranes[interfaces[k].sourceDomain]->z.col(0);
                        Vtsource(k) = h_1s1%V1.col(0);
                        Vnsource(k) = h_2s1%V1.col(0) + h_2s2%V2.col(0);
                        break;
                    }
                    case 1: // East
                    {
                        size_t nx = membranes[interfaces[k].sourceDomain]->nx;
                        arma::rowvec h_1s1 = sqrt(e11.row(nx-1));
                        arma::rowvec h_1s2 = e12.row(nx-1)/sqrt(e11.row(nx-1));
                        arma::rowvec h_2s2 = 1/sqrt(e_22.row(nx-1));
                        Zsource(k) = membranes[interfaces[k].sourceDomain]->z.row(nx-1).t();
                        Vnsource(k) = (h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t();
                        Vtsource(k) = (h_2s2%V2.row(nx-1)).t();
                        break;
                    }
                    case 2: // North
                    {
                        size_t ny = membranes[interfaces[k].sourceDomain]->ny;
                        arma::vec h_1s1 = 1/sqrt(e_11.col(ny-1));
                        arma::vec h_2s1 = e12.col(ny-1)/sqrt(e22.col(ny-1));
                        arma::vec h_2s2 = sqrt(e22.col(ny-1));
                        Zsource(k) = membranes[interfaces[k].sourceDomain]->z.col(ny-1);
                        Vtsource(k) = h_1s1%V1.col(ny-1);
                        Vnsource(k) = h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1);
                        break;
                    }
                    case 3: // West
                    {
                        arma::rowvec h_1s1 = sqrt(e11.row(0));
                        arma::rowvec h_1s2 = e12.row(0)/sqrt(e11.row(0));
                        arma::rowvec h_2s2 = 1/sqrt(e_22.row(0));
                        Zsource(k) = membranes[interfaces[k].sourceDomain]->z.row(0).t();
                        Vnsource(k) = (h_1s1%V1.row(0) + h_1s2%V2.row(0)).t();
                        Vtsource(k) = (h_2s2%V2.row(0)).t();
                        break;
                    }
                }
            }
            converged = true;
            for (size_t k = 0; k < interfaces.size(); k++)
            {
                printf("Interface %lu\n", k+1);
                double res1 = arma::norm(Vntarget(k) - Vnsource(k))/arma::norm(Vntarget(k) + Vnsource(k));
                double res2 = arma::norm(Vttarget(k) - Vtsource(k))/arma::norm(Vttarget(k) + Vtsource(k));
                double res3 = arma::norm(Ztarget(k)  - Zsource(k)) /arma::norm(Ztarget(k)  + Zsource(k));
                printf("Residual normal  %4.2e\n", res1);
                printf("Residual tangent %4.2e\n", res2);
                printf("Residual z       %4.2e\n", res3);
                if (res1 > residualTarget || res2 > residualTarget || res3 > residualTarget)
                    converged = false;
            }
            if (count > iterations)
                converged = true;
        } while (converged == false);
    }
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->iter = iter_old(k);
}