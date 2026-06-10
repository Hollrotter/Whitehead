#include "Structure.hpp"

void Structure::linear()
{
    analysis = Analysis::linear;
    bool converged = false;
    size_t count = 1;
    arma::field<arma::vec> Ztarget(interfaces.size()), Zsource(interfaces.size());
    for (const Interface& interface:interfaces)
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
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.targetCurve)
                {
                    case 0: // South
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
                        break;
                    case 1: // East
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 2: // North
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource, 1);
                        break;
                    case 3: // West
                        membraneTarget->boundary(Field::z, targetDirection, BC::Robin, interface.lambdaSource,-1);
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
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 1: // East
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 2: // North
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
            case 3: // West
            {
                switch (interface.sourceCurve)
                {
                    case 0: // South
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                    case 1: // East
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 2: // North
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget, 1);
                        break;
                    case 3: // West
                        membraneSource->boundary(Field::z, sourceDirection, BC::Robin, interface.lambdaTarget,-1);
                        break;
                }
                break;
            }
        }
    }
    #pragma omp parallel for
    for (auto& membrane:membranes)
        membrane->solve_S();
    do
    {
        std::cout << "Iteration " << count << '/' << iterations << std::endl;
        for (const Interface& interface:interfaces)
        {
            Membrane *membraneSource = membranes[interface.sourceDomain];
            Membrane *membraneTarget = membranes[interface.targetDomain];
            size_t nx = membraneSource->nx;
            size_t ny = membraneSource->ny;
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
                            sourceZ = reverse(sourceZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i) = sourceZ(i);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1) = sourceZ(nx-1);
                            break;
                        case 1: // East
                            sourceZ = reverse(sourceZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(membraneTarget->nx-1+j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(membraneTarget->nx-1) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                                membraneTarget->b(membraneTarget->nx-1+(ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i+(membraneTarget->ny-1)*nx) = sourceZ(i);
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((membraneTarget->ny-1)*nx) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1+(membraneTarget->ny-1)*nx) = sourceZ(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneSource->h_1s1_east;
                    arma::rowvec h_1s2 = membraneSource->h_1s2_east;
                    arma::vec sourceZ = (interface.lambdaSource*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            sourceZ = reverse(sourceZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i) = sourceZ(i);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1) = sourceZ(nx-1);
                            break;
                        case 1: // East
                            sourceZ = reverse(sourceZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(membraneTarget->nx-1+j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(membraneTarget->nx-1) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                                membraneTarget->b(membraneTarget->nx-1+(ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i+(membraneTarget->ny-1)*nx) = sourceZ(i);
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((membraneTarget->ny-1)*nx) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1+(membraneTarget->ny-1)*nx) = sourceZ(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneSource->h_2s1_north;
                    arma::vec h_2s2 = membraneSource->h_2s2_north;
                    arma::vec sourceZ = interface.lambdaSource*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.targetCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i) = sourceZ(i);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1) = sourceZ(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(membraneTarget->nx-1+j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(membraneTarget->nx-1) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                                membraneTarget->b(membraneTarget->nx-1+(ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                        case 2: // North
                            sourceZ = reverse(sourceZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i+(membraneTarget->ny-1)*nx) = sourceZ(i);
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((membraneTarget->ny-1)*nx) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1+(membraneTarget->ny-1)*nx) = sourceZ(nx-1);
                            break;
                        case 3: // West
                            sourceZ = reverse(sourceZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((ny-1)*membraneTarget->nx) = sourceZ(ny-1);
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
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i) = sourceZ(i);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[0]->curveType == CurveType::Boundary  && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1) = sourceZ(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(membraneTarget->nx-1+j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(membraneTarget->nx-1) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Boundary)
                                membraneTarget->b(membraneTarget->nx-1+(ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                        case 2: // North
                            sourceZ = reverse(sourceZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneTarget->b(i+(membraneTarget->ny-1)*nx) = sourceZ(i);
                            if (membraneTarget->chi[2]->curveType == CurveType::Boundary  && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((membraneTarget->ny-1)*nx) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[1]->curveType == CurveType::Interface)
                                membraneTarget->b(nx-1+(membraneTarget->ny-1)*nx) = sourceZ(nx-1);
                            break;
                        case 3: // West
                            sourceZ = reverse(sourceZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneTarget->b(j*membraneTarget->nx) = sourceZ(j);
                            if (membraneTarget->chi[0]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Boundary)
                                membraneTarget->b(0) = sourceZ(0);
                            if (membraneTarget->chi[2]->curveType == CurveType::Interface && membraneTarget->chi[3]->curveType == CurveType::Interface)
                                membraneTarget->b((ny-1)*membraneTarget->nx) = sourceZ(ny-1);
                            break;
                    }
                    break;
                }
            }
            nx = membraneTarget->nx;
            ny = membraneTarget->ny;
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
                            targetZ = reverse(targetZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i) = targetZ(i);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1) = targetZ(nx-1);
                            break;
                        case 1: // East
                            targetZ = reverse(targetZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(membraneSource->nx-1+j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(membraneSource->nx-1) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                                membraneSource->b(membraneSource->nx-1+(ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i+(membraneSource->ny-1)*nx) = targetZ(i);
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((membraneSource->ny-1)*nx) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1+(membraneSource->ny-1)*nx) = targetZ(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                    }
                    break;
                }
                case 1: // East
                {
                    arma::rowvec dZd1 = D1.row(nx-1)*Z;
                    arma::rowvec dZd2 = Z.row(nx-1)*D2.t();
                    arma::rowvec h_1s1 = membraneTarget->h_1s1_east;
                    arma::rowvec h_1s2 = membraneTarget->h_1s2_east;
                    arma::vec targetZ = (interface.lambdaTarget*Z.row(nx-1) - (h_1s1%dZd1 + h_1s2%dZd2)).t();
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            targetZ = reverse(targetZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i) = targetZ(i);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1) = targetZ(nx-1);
                            break;
                        case 1: // East
                            targetZ = reverse(targetZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(membraneSource->nx-1+j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(membraneSource->nx-1) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                                membraneSource->b(membraneSource->nx-1+(ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                        case 2: // North
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i+(membraneSource->ny-1)*nx) = targetZ(i);
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((membraneSource->ny-1)*nx) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1+(membraneSource->ny-1)*nx) = targetZ(nx-1);
                            break;
                        case 3: // West
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                    }
                    break;
                }
                case 2: // North
                {
                    arma::vec dZd1 = D1*Z.col(ny-1);
                    arma::vec dZd2 = Z*D2.row(ny-1).t();
                    arma::vec h_2s1 = membraneTarget->h_2s1_north;
                    arma::vec h_2s2 = membraneTarget->h_2s2_north;
                    arma::vec targetZ = interface.lambdaTarget*Z.col(ny-1) - (h_2s1%dZd1 + h_2s2%dZd2);
                    switch (interface.sourceCurve)
                    {
                        case 0: // South
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i) = targetZ(i);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1) = targetZ(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(membraneSource->nx-1+j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(membraneSource->nx-1) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                                membraneSource->b(membraneSource->nx-1+(ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                        case 2: // North
                            targetZ = reverse(targetZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i+(membraneSource->ny-1)*nx) = targetZ(i);
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((membraneSource->ny-1)*nx) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1+(membraneSource->ny-1)*nx) = targetZ(nx-1);
                            break;
                        case 3: // West
                            targetZ = reverse(targetZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((ny-1)*membraneSource->nx) = targetZ(ny-1);
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
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i) = targetZ(i);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[0]->curveType == CurveType::Boundary  && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1) = targetZ(nx-1);
                            break;
                        case 1: // East
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(membraneSource->nx-1+j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(membraneSource->nx-1) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Boundary)
                                membraneSource->b(membraneSource->nx-1+(ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                        case 2: // North
                            targetZ = reverse(targetZ);
                            for (size_t i = 1; i < nx-1; i++)
                                membraneSource->b(i+(membraneSource->ny-1)*nx) = targetZ(i);
                            if (membraneSource->chi[2]->curveType == CurveType::Boundary  && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((membraneSource->ny-1)*nx) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[1]->curveType == CurveType::Interface)
                                membraneSource->b(nx-1+(membraneSource->ny-1)*nx) = targetZ(nx-1);
                            break;
                        case 3: // West
                            targetZ = reverse(targetZ);
                            for (size_t j = 1; j < ny-1; j++)
                                membraneSource->b(j*membraneSource->nx) = targetZ(j);
                            if (membraneSource->chi[0]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Boundary)
                                membraneSource->b(0) = targetZ(0);
                            if (membraneSource->chi[2]->curveType == CurveType::Interface && membraneSource->chi[3]->curveType == CurveType::Interface)
                                membraneSource->b((ny-1)*membraneSource->nx) = targetZ(ny-1);
                            break;
                    }
                    break;
                }
            }
        }
        #pragma omp parallel for
        for (auto& membrane:membranes)
            membrane->solve_b();

        for (size_t k = 0; k < interfaces.size(); k++)
        {
            Interface interface = interfaces[k];
            const Membrane *membraneSource = membranes[interface.sourceDomain];
            arma::mat Z  = membraneSource->z;
            arma::mat D1 = membraneSource->D1;
            arma::mat D2 = membraneSource->D2;
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
            Z  = membraneTarget->z;
            D1 = membraneTarget->D1;
            D2 = membraneTarget->D2;
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
            std::cout << "Interface " << k+1 << "\n";
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