#include "Structure.hpp"

void Structure::planeStrain()
{
    bool converged = false;
    int count = 1;
    arma::field<arma::vec> Vntarget(interfaces.size()), Vnsource(interfaces.size());
    arma::field<arma::vec> Vttarget(interfaces.size()), Vtsource(interfaces.size());
    for (Interface& interface:interfaces)
    {
        Direction targetDirection = static_cast<Direction>(interface.targetCurve);
        Membrane *membraneTarget = membranes[interface.targetDomain];
        switch (interface.sourceCurve)
        {
            case 0: // South
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource,-1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::v1, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        membraneTarget->boundary(Field::v2, targetDirection, BC::Robin,-interface.lambdaSource, 1);
                        break;
                }
                break;
            }
        }
        Direction sourceDirection = static_cast<Direction>(interface.sourceCurve);
        Membrane *membraneSource = membranes[interface.sourceDomain];
        switch (interface.targetCurve)
        {
            case 0: // South
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget,-1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::v1, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        membraneSource->boundary(Field::v2, sourceDirection, BC::Robin,-interface.lambdaTarget, 1);
                        break;
                }
                break;
            }
        }
    }
    #pragma omp parallel for
    for (auto& membrane:membranes)
        membrane->planeStrainSolve();
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << '\n';
        for (Interface& interface:interfaces)
        {
            Membrane *membraneSource = membranes[interface.sourceDomain];
            Membrane *membraneTarget = membranes[interface.targetDomain];
            size_t    nx    = membraneSource->nx;
            size_t    ny    = membraneSource->ny;
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
                    arma::vec h_1s1 = membraneSource->h_1s1_south;
                    arma::vec h_2s1 = membraneSource->h_2s1_south;
                    arma::vec h_2s2 = membraneSource->h_2s2_south;
                    arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(0)
                                       + interface.c1Source%v1__1.col(0) + interface.c2Source%v1__2.col(0);
                    arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                       + interface.c11Source%v1__1.col(0) + interface.c22Source%v2__2.col(0) + interface.c12Source%(v1__2.col(0) + v2__1.col(0));
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i)                     = sourceV1(i);
                                membraneTarget->bv(i+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 1: // East
                            std::swap(sourceV1, sourceV2);
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(membraneTarget->nx-1)                     = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nx-1+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx)                     = sourceV1(i);
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((membraneTarget->ny-1)*nx)                     = sourceV1(0);
                                membraneTarget->bv((membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 3: // West
                            std::swap(sourceV1, sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv((ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                    arma::rowvec h_2s2 = membraneSource->h_2s2_east;
                    arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                       - interface.c11Source%v1__1.row(nx-1).t() - interface.c22Source%v2__2.row(nx-1).t() - interface.c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                       - interface.c1Source%v2__1.row(nx-1).t() - interface.c2Source%v2__2.row(nx-1).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            std::swap(sourceV1, sourceV2);
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i)                     = sourceV1(i);
                                membraneTarget->bv(i+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 1: // East
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(membraneTarget->nx-1)                     = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nx-1+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                        case 2: // North
                            std::swap(sourceV1, sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx)                     = sourceV1(i);
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((membraneTarget->ny-1)*nx)                     = sourceV1(0);
                                membraneTarget->bv((membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv((ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec h_1s1 = membraneSource->h_1s1_north;
                    arma::vec h_2s1 = membraneSource->h_2s1_north;
                    arma::vec h_2s2 = membraneSource->h_2s2_north;
                    arma::vec sourceV1 = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                       - interface.c1Source%v1__1.col(ny-1) - interface.c2Source%v1__2.col(ny-1);
                    arma::vec sourceV2 = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                       - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i)                     = sourceV1(i);
                                membraneTarget->bv(i+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 1: // East
                            std::swap(sourceV1, sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(membraneTarget->nx-1)                     = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nx-1+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                        case 2: // North
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx)                     = sourceV1(i);
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((membraneTarget->ny-1)*nx)                     = sourceV1(0);
                                membraneTarget->bv((membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 3: // West
                            std::swap(sourceV1, sourceV2);
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv((ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membraneSource->h_1s1_west;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_west;
                    arma::rowvec h_2s2 = membraneSource->h_2s2_west;
                    arma::vec sourceV1 = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                       + interface.c11Source%v1__1.row(0).t() + interface.c22Source%v2__2.row(0).t() + interface.c12Source%(v1__2.row(0) + v2__1.row(0)).t();
                    arma::vec sourceV2 = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                       + interface.c1Source%v2__1.row(0).t() + interface.c2Source%v2__2.row(0).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            std::swap(sourceV1, sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i)                     = sourceV1(i);
                                membraneTarget->bv(i+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(membraneTarget->nx-1+j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(membraneTarget->nx-1)                     = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nx-1+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv(membraneTarget->nx-1+(ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                        case 2: // North
                            std::swap(sourceV1, sourceV2);
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx)                     = sourceV1(i);
                                membraneTarget->bv(i+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(i);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((membraneTarget->ny-1)*nx)                     = sourceV1(0);
                                membraneTarget->bv((membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx)                     = sourceV1(nx-1);
                                membraneTarget->bv(nx-1+(membraneTarget->ny-1)*nx+membraneTarget->nxy) = sourceV2(nx-1);
                            }
                            break;
                        case 3: // West
                            sourceV1 = reverse(sourceV1);
                            sourceV2 = reverse(sourceV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneTarget->bv(j*membraneTarget->nx)                     = sourceV1(j);
                                membraneTarget->bv(j*membraneTarget->nx+membraneTarget->nxy) = sourceV2(j);
                            }
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneTarget->bv(0)                   = sourceV1(0);
                                membraneTarget->bv(membraneTarget->nxy) = sourceV2(0);
                            }
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneTarget->bv((ny-1)*membraneTarget->nx)                     = sourceV1(ny-1);
                                membraneTarget->bv((ny-1)*membraneTarget->nx+membraneTarget->nxy) = sourceV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
            }
            nx    = membraneTarget->nx;
            ny    = membraneTarget->ny;
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
                    arma::vec h_1s1 = membraneTarget->h_1s1_south;
                    arma::vec h_2s1 = membraneTarget->h_2s1_south;
                    arma::vec h_2s2 = membraneTarget->h_2s2_south;
                    arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(0)
                                       + interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0);
                    arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                       + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i)                     = targetV1(i);
                                membraneSource->bv(i+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 1: // East
                            std::swap(targetV1, targetV2);
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(membraneSource->nx-1)                     = targetV1(0);
                                membraneSource->bv(membraneSource->nx-1+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i+(membraneSource->ny-1)*nx)                     = targetV1(i);
                                membraneSource->bv(i+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((membraneSource->ny-1)*nx)                     = targetV1(0);
                                membraneSource->bv((membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 3: // West
                            std::swap(targetV1, targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv((ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_east;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_east;
                    arma::rowvec h_2s2 = membraneTarget->h_2s2_east;
                    arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                       - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(nx-1)).t()
                                       - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            std::swap(targetV1, targetV2);
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i)                     = targetV1(i);
                                membraneSource->bv(i+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 1: // East
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(membraneSource->nx-1)                     = targetV1(0);
                                membraneSource->bv(membraneSource->nx-1+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                        case 2: // North
                            std::swap(targetV1, targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i+(membraneSource->ny-1)*nx)                     = targetV1(i);
                                membraneSource->bv(i+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((membraneSource->ny-1)*nx)                     = targetV1(0);
                                membraneSource->bv((membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv((ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec h_1s1 = membraneTarget->h_1s1_north;
                    arma::vec h_2s1 = membraneTarget->h_2s1_north;
                    arma::vec h_2s2 = membraneTarget->h_2s2_north;
                    arma::vec targetV1 = interface.lambdaTarget*h_1s1%V1.col(ny-1)
                                       - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                    arma::vec targetV2 = interface.lambdaTarget*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                       - interface.c11Target%v1__1.col(ny-1) - interface.c22Target%v2__2.col(ny-1) - interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i)                     = targetV1(i);
                                membraneSource->bv(i+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 1: // East
                            std::swap(targetV1, targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(membraneSource->nx-1)                     = targetV1(0);
                                membraneSource->bv(membraneSource->nx-1+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                        case 2: // North
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i+(membraneSource->ny-1)*nx)                     = targetV1(i);
                                membraneSource->bv(i+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((membraneSource->ny-1)*nx)                     = targetV1(0);
                                membraneSource->bv((membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 3: // West
                            std::swap(targetV1, targetV2);
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv((ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_west;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_west;
                    arma::rowvec h_2s2 = membraneTarget->h_2s2_west;
                    arma::vec targetV1 = interface.lambdaTarget*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                       + interface.c11Target%v1__1.row(0).t() + interface.c22Target%v2__2.row(0).t() + interface.c12Target%(v1__2.row(0) + v2__1.row(0)).t();
                    arma::vec targetV2 = interface.lambdaTarget*(h_2s2%V2.row(0)).t()
                                       + interface.c1Target%v2__1.row(0).t() + interface.c2Target%v2__2.row(0).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            std::swap(targetV1, targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i)                     = targetV1(i);
                                membraneSource->bv(i+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(membraneSource->nx-1+j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(membraneSource->nx-1)                     = targetV1(0);
                                membraneSource->bv(membraneSource->nx-1+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv(membraneSource->nx-1+(ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                        case 2: // North
                            std::swap(targetV1, targetV2);
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t i = 1; i < nx-1; i++)
                            {
                                membraneSource->bv(i+(membraneSource->ny-1)*nx)                     = targetV1(i);
                                membraneSource->bv(i+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(i);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((membraneSource->ny-1)*nx)                     = targetV1(0);
                                membraneSource->bv((membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx)                     = targetV1(nx-1);
                                membraneSource->bv(nx-1+(membraneSource->ny-1)*nx+membraneSource->nxy) = targetV2(nx-1);
                            }
                            break;
                        case 3: // West
                            targetV1 = reverse(targetV1);
                            targetV2 = reverse(targetV2);
                            for (size_t j = 1; j < ny-1; j++)
                            {
                                membraneSource->bv(j*membraneSource->nx)                     = targetV1(j);
                                membraneSource->bv(j*membraneSource->nx+membraneSource->nxy) = targetV2(j);
                            }
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                            {
                                membraneSource->bv(0)                   = targetV1(0);
                                membraneSource->bv(membraneSource->nxy) = targetV2(0);
                            }
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                            {
                                membraneSource->bv((ny-1)*membraneSource->nx)                     = targetV1(ny-1);
                                membraneSource->bv((ny-1)*membraneSource->nx+membraneSource->nxy) = targetV2(ny-1);
                            }
                            break;
                    }
                    break;
                }
            }
        }
        #pragma omp parallel for
        for (auto& membrane:membranes)
            membrane->planeStrainEval();
        for (size_t k = 0; k < interfaces.size(); k++)
        {
            Interface interface = interfaces[k];
            Membrane *membraneSource = membranes[interface.sourceDomain];
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
                    arma::vec h_1s1 = membraneSource->h_1s1_south;
                    arma::vec h_2s1 = membraneSource->h_2s1_south;
                    arma::vec h_2s2 = membraneSource->h_2s2_south;
                    Vtsource(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                + interface.c1Source%v1__1.col(0) + interface.c2Source%v1__2.col(0);
                    Vnsource(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                + interface.c11Source%v1__1.col(0) + interface.c22Source%v2__2.col(0) + interface.c12Source%(v1__2.col(0) + v2__1.col(0));
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneSource->nx;
                    arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                    arma::rowvec h_2s2 = membraneSource->h_2s2_east;
                    Vnsource(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                - interface.c11Source%v1__1.row(nx-1).t() - interface.c22Source%v2__2.row(nx-1).t() - interface.c12Source%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    Vtsource(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                - interface.c1Source%v2__1.row(nx-1).t() - interface.c2Source%v2__2.row(nx-1).t();
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneSource->ny;
                    arma::vec h_1s1 = membraneSource->h_1s1_north;
                    arma::vec h_2s1 = membraneSource->h_2s1_north;
                    arma::vec h_2s2 = membraneSource->h_2s2_north;
                    Vtsource(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                - interface.c1Source%v1__1.col(ny-1) - interface.c2Source%v1__2.col(ny-1);
                    Vnsource(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                - interface.c11Source%v1__1.col(ny-1) - interface.c22Source%v2__2.col(ny-1) - interface.c12Source%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membraneSource->h_1s1_west;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_west;
                    arma::rowvec h_2s2 = membraneSource->h_2s2_west;
                    Vnsource(k) = interface.lambdaSource*(h_1s1%V1.row(0) + h_1s2%V2.row(0)).t()
                                + interface.c11Source%v1__1.row(0).t() + interface.c22Source%v2__2.row(0).t() + interface.c12Source%(v1__2.row(0) + v2__1.row(0)).t();
                    Vtsource(k) = interface.lambdaSource*(h_2s2%V2.row(0)).t()
                                + interface.c1Source%v2__1.row(0).t() + interface.c2Source%v2__2.row(0).t();
                    break;
                }
            }
            Membrane *membraneTarget = membranes[interface.targetDomain];
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
                    arma::vec h_1s1 = membraneTarget->h_1s1_south;
                    arma::vec h_2s1 = membraneTarget->h_2s1_south;
                    arma::vec h_2s2 = membraneTarget->h_2s2_south;
                    if (interface.sourceCurve == 0) // South
                    {
                        Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(0)
                                    + (interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0));
                        Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    + (interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0)));
                    }
                    else if (interface.sourceCurve == 1) // East
                    {
                        Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(0)
                                    + (interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0));
                        Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    - (interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0)));
                    }
                    else if (interface.sourceCurve == 2) // North
                    {
                        Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                    - (interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0));
                        Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    - (interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0)));
                    }
                    else // West
                    {
                        Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(0)
                                    - (interface.c1Target%v1__1.col(0) + interface.c2Target%v1__2.col(0));
                        Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(0) + h_2s2%V2.col(0))
                                    + interface.c11Target%v1__1.col(0) + interface.c22Target%v2__2.col(0) + interface.c12Target%(v1__2.col(0) + v2__1.col(0));
                    }
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                    {
                        Vttarget(k) = reverse(Vttarget(k));
                        Vntarget(k) = reverse(Vntarget(k));
                    }
                    break;
                }
                case 1: // East
                {
                    size_t nx = membraneTarget->nx;
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_east;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_east;
                    arma::rowvec h_2s2 = membraneTarget->h_2s2_east;
                    if (interface.sourceCurve == 0) // South
                    {
                        Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    + interface.c11Target%v1__1.row(nx-1).t() + interface.c22Target%v2__2.row(nx-1).t() + interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    }
                    else if (interface.sourceCurve == 1) // East
                    {
                        Vttarget(k) =-interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    - interface.c1Target%v2__1.row(nx-1).t() - interface.c2Target%v2__2.row(nx-1).t();
                        Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    }
                    else if (interface.sourceCurve == 2) // North
                    {
                        Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    + interface.c1Target%v2__1.row(nx-1).t() + interface.c2Target%v2__2.row(nx-1).t();
                        Vntarget(k) =-interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    - interface.c11Target%v1__1.row(nx-1).t() - interface.c22Target%v2__2.row(nx-1).t() - interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    }
                    else // West
                    {
                        Vttarget(k) = interface.lambdaSource*(h_2s2%V2.row(nx-1)).t()
                                    + interface.c1Target%v2__1.row(nx-1).t() + interface.c2Target%v2__2.row(nx-1).t();
                        Vntarget(k) = interface.lambdaSource*(h_1s1%V1.row(nx-1) + h_1s2%V2.row(nx-1)).t()
                                    + interface.c11Target%v1__1.row(nx-1).t() + interface.c22Target%v2__2.row(nx-1).t() + interface.c12Target%(v1__2.row(nx-1) + v2__1.row(nx-1)).t();
                    }
                    if (interface.sourceCurve == 0 || interface.sourceCurve == 1)
                    {
                        Vttarget(k) = reverse(Vttarget(k));
                        Vntarget(k) = reverse(Vntarget(k));
                    }
                    break;
                }
                case 2: // North
                {
                    size_t ny = membraneTarget->ny;
                    arma::vec h_1s1 = membraneTarget->h_1s1_north;
                    arma::vec h_2s1 = membraneTarget->h_2s1_north;
                    arma::vec h_2s2 = membraneTarget->h_2s2_north;
                    if (interface.sourceCurve == 0) // South
                    {
                        Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                    + interface.c1Target%v1__1.col(ny-1) + interface.c2Target%v1__2.col(ny-1);
                        Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    + interface.c11Target%v1__1.col(ny-1) + interface.c22Target%v2__2.col(ny-1) + interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    }
                    else if (interface.sourceCurve == 1) // East
                    {
                        Vttarget(k) = interface.lambdaSource*h_1s1%V1.col(ny-1)
                                    + interface.c1Target%v1__1.col(ny-1) + interface.c2Target%v1__2.col(ny-1);
                        Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    - interface.c11Target%v1__1.col(ny-1) - interface.c22Target%v2__2.col(ny-1) - interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    }
                    else if (interface.sourceCurve == 2) // North
                    {
                        Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(ny-1)
                                    - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                        Vntarget(k) =-interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    - interface.c11Target%v1__1.col(ny-1) - interface.c22Target%v2__2.col(ny-1) - interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    }
                    else // West
                    {
                        Vttarget(k) =-interface.lambdaSource*h_1s1%V1.col(ny-1)
                                    - interface.c1Target%v1__1.col(ny-1) - interface.c2Target%v1__2.col(ny-1);
                        Vntarget(k) = interface.lambdaSource*(h_2s1%V1.col(ny-1) + h_2s2%V2.col(ny-1))
                                    + interface.c11Target%v1__1.col(ny-1) + interface.c22Target%v2__2.col(ny-1) + interface.c12Target%(v1__2.col(ny-1) + v2__1.col(ny-1));
                    }
                    if (interface.sourceCurve == 2 || interface.sourceCurve == 3)
                    {
                        Vttarget(k) = reverse(Vttarget(k));
                        Vntarget(k) = reverse(Vntarget(k));
                    }
                    break;
                }
                case 3: // West
                {
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_west;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_west;
                    arma::rowvec h_2s2 = membraneTarget->h_2s2_west;
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
            printf("Residual normal  %4.2e\n", res1);
            printf("Residual tangent %4.2e\n", res2);
            if (res1 > residualTarget || res2 > residualTarget)
                converged = false;
        }
        count++;
        if (count > iterations)
            converged = true;
    } while (converged == false);
}