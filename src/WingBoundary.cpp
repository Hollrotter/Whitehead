#include "Wing.hpp"

template <class C> void Wing::boundary(const Direction dir, const BC bc, const C val)
{
    switch (dir)
    {
        case Direction::N:
            mu.northBC = bc;
            mu.north   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.north*val;
            break;
        case Direction::S:
            mu.southBC = bc;
            mu.south   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.south*val;
            break;
        case Direction::W:
            mu.westBC = bc;
            mu.west   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.west*val;
            break;
        case Direction::E:
            mu.eastBC = bc;
            mu.east   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.east*val;
            break;
    }
}

template <class C> void Wing::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const C val)
{
    if (dir == chi[0])
    {
        mu.southBC = bc;
        mu.south   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.south*val;
    }
    else if (dir == chi[1])
    {
        mu.eastBC = bc;
        mu.east   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.east*val;
    }
    else if (dir == chi[2])
    {
        mu.northBC = bc;
        mu.north   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.north*val;
    }
    else if (dir == chi[3])
    {
        mu.westBC = bc;
        mu.west   = typeid(val) == typeid(arma::vec) ? arma::vec(val) : mu.west*val;
    }
    else
    {
        std::println("Not a curve of the wing!");
        exit(EXIT_FAILURE);
    }
}

template <class C> void Wing::boundary(const Direction dir, const BC bc, const double _r1, const double _r2, const C val)
{
    switch (dir)
    {
        case Direction::N:
            mu.r1North = _r1;
            mu.r2North = _r2;
            break;
        case Direction::S:
            mu.r1South = _r1;
            mu.r2South = _r2;
            break;
        case Direction::W:
            mu.r1West = _r1;
            mu.r2West = _r2;
            break;
        case Direction::E:
            mu.r1East = _r1;
            mu.r2East = _r2;
            break;
    }
    boundary(dir, bc, val);
}

template <class C> void Wing::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2, const C val)
{
    if (dir == chi[0])
    {
        mu.r1South = _r1;
        mu.r2South = _r2;
    }
    else if (dir == chi[1])
    {
        mu.r1East = _r1;
        mu.r2East = _r2;
    }
    else if (dir == chi[2])
    {
        mu.r1North = _r1;
        mu.r2North = _r2;
    }
    else if (dir == chi[3])
    {
        mu.r1West = _r1;
        mu.r2West = _r2;
    }
    else
    {
        std::println("Not a curve of the wing!");
        exit(EXIT_FAILURE);
    }
    boundary(dir, bc, val);
}

void Wing::muBoundarySouth(const size_t i)
{
    switch (mu.southBC)
    {
        case BC::Dirichlet:
        {
            double w1 = phi1->weightFunction(xi_1(i));
            double w2 = phi2->weightFunction(-1);
            for (size_t q = 0; q < ny; q++)
            {
                double t2 = phi2->left(q);
                for (size_t p = 0; p < nx; p++)
                    A(i, p+q*nx) = w1*PHI1(i, p) * w2*t2;
            }
            b(i) = mu.south(i);
            break;
        }
        case BC::Neumann:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(-1);
            double dw2 = phi2->weightFunctionDerivative(-1);
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->left(q);
                double dt2 = phi2->leftDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(i, p+q*nx) = h_2s1_south(i)*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                 + h_2s2_south(i)*w1*PHI1(i, p)*(w2*dt2 + dw2*t2);
            }
            b(i) = mu.south(i);
            break;
        }
        case BC::Robin:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(-1);
            double dw2 = phi2->weightFunctionDerivative(-1);
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->left(q);
                double dt2 = phi2->leftDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(i, p+q*nx) = mu.r1South*w1*PHI1(i, p) * w2*t2
                                 + mu.r2South*(h_2s1_south(i)*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                             + h_2s2_south(i)*w1*PHI1(i, p)*(w2*dt2 + dw2*t2));
            }
            b(i) = mu.south(i);
            break;
        }
        case BC::Kutta:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(-1);
            double dw2 = phi2->weightFunctionDerivative(-1);
            auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), -1, chi);
            double detJ = dxdx1*dydx2 - dxdx2*dydx1;
            double J11_inv = dydx2/detJ;
            double J21_inv =-dydx1/detJ;
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->left(q);
                double dt2 = phi2->leftDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(i, p+q*nx) = J11_inv*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                 + J21_inv*w1*PHI1(i, p)*(w2*dt2 + dw2*t2);
            }
            b(i) = mu.south(i);
            break;
        }
        case BC::None:
            break;
    }
}

void Wing::muBoundaryNorth(const size_t i)
{
    size_t k = i+(ny-1)*nx;
    switch (mu.northBC)
    {
        case BC::Dirichlet:
        {
            double w1 = phi1->weightFunction(xi_1(i));
            double w2 = phi2->weightFunction(1);
            for (size_t q = 0; q < ny; q++)
            {
                double t2 = phi2->right(q);
                for (size_t p = 0; p < nx; p++) 
                    A(k, p+q*nx) = w1*PHI1(i, p) * w2*t2;
            }
            b(k) = mu.north(i);
            break;
        }
        case BC::Neumann:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(1);
            double dw2 = phi2->weightFunctionDerivative(1);
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->right(q);
                double dt2 = phi2->rightDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = h_2s1_north(i)*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                 + h_2s2_north(i)*w1*PHI1(i, p)*(w2*dt2 + dw2*t2);
            }
            b(k) = mu.north(i);
            break;
        }
        case BC::Robin:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(1);
            double dw2 = phi2->weightFunctionDerivative(1);
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->right(q);
                double dt2 = phi2->rightDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = mu.r1North*w1*PHI1(i, p) * w2*t2
                                 + mu.r2North*(h_2s1_north(i)*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                             + h_2s2_north(i)*w1*PHI1(i, p)*(w2*dt2 + dw2*t2));
            }
            b(k) = mu.north(i);
            break;
        }
        case BC::Kutta:
        {
            double  w1 = phi1->weightFunction(xi_1(i));
            double dw1 = phi1->weightFunctionDerivative(xi_1(i));
            double  w2 = phi2->weightFunction(1);
            double dw2 = phi2->weightFunctionDerivative(1);
            auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), 1, chi);
            double detJ = dxdx1*dydx2 - dxdx2*dydx1;
            double J11_inv = dydx2/detJ;
            double J21_inv =-dydx1/detJ;
            for (size_t q = 0; q < ny; q++)
            {
                double  t2 = phi2->right(q);
                double dt2 = phi2->rightDerivative(q);
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = J11_inv*(w1*dPHI1(i, p) + dw1*PHI1(i, p))*w2*t2
                                 + J21_inv*w1*PHI1(i, p)*(w2*dt2 + dw2*t2);
            }
            b(k) = mu.north(i);
            break;
        }
        case BC::None:
            break;
    }
}

void Wing::muBoundaryWest(const size_t j)
{
    size_t k = j*nx;
    switch (mu.westBC)
    {
        case BC::Dirichlet:
        {
            double w1 = phi1->weightFunction(-1);
            double w2 = phi2->weightFunction(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double t1 = phi1->left(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = w1*t1 * w2*PHI2(j, q);
            }
            b(k) = mu.west(j);
            break;
        }
        case BC::Neumann:
        {
            double  w1 = phi1->weightFunction(-1);
            double dw1 = phi1->weightFunctionDerivative(-1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->left(p);
                double dt1 = phi1->leftDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = h_1s1_west(j)*(w1*dt1 + dw1*t1)*w2*PHI2(j, q)
                                 + h_1s2_west(j)*w1*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q));
            }
            b(k) = mu.west(j);
            break;
        }
        case BC::Robin:
        {
            double  w1 = phi1->weightFunction(-1);
            double dw1 = phi1->weightFunctionDerivative(-1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->left(p);
                double dt1 = phi1->leftDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = mu.r1West*w1*t1 * w2*PHI2(j, q)
                                 + mu.r2West*(h_1s1_west(j)*(w1*dt1 + dw1*t1)*w2*PHI2(j, q)
                                            + h_1s2_west(j)*w1*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q)));
            }
            b(k) = mu.west(j);
            break;
        }
        case BC::Kutta:
        {
            double  w1 = phi1->weightFunction(-1);
            double dw1 = phi1->weightFunctionDerivative(-1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(-1, xi_2(j), chi);
            double detJ = dxdx1*dydx2 - dxdx2*dydx1;
            double J11_inv = dydx2/detJ;
            double J21_inv =-dydx1/detJ;
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->left(p);
                double dt1 = phi1->leftDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = J11_inv*(w1*dt1 + dw1*t1)*w2*PHI2(j, q)
                                 + J21_inv*w1*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q));
            }
            b(k) = mu.west(j);
            break;
        }
        case BC::None:
            break;
    }
}

void Wing::muBoundaryEast(const size_t j)
{
    size_t k = nx-1+j*nx;
    switch (mu.eastBC)
    {
        case BC::Dirichlet:
        {
            double w1 = phi1->weightFunction(1);
            double w2 = phi2->weightFunction(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double t1 = phi1->right(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = w1*t1 * w2*PHI2(j, q);
            }
            b(k) = mu.east(j);
            break;
        }
        case BC::Neumann:
        {
            double  w1 = phi1->weightFunction(1);
            double dw1 = phi1->weightFunctionDerivative(1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->right(p);
                double dt1 = phi1->rightDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = h_1s1_east(j)*(w1*dt1 + dw1*t1)*w2*PHI2(j, q)
                                 + h_1s2_east(j)*w1*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q));
            }
            b(k) = mu.east(j);
            break;
        }
        case BC::Robin:
        {
            double  w1 = phi1->weightFunction(1);
            double dw1 = phi1->weightFunctionDerivative(1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->right(p);
                double dt1 = phi1->rightDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = mu.r1East*w1*t1 * w2*PHI2(j, q)
                                 + mu.r2East*(h_1s1_east(j)*(w1*dt1 + dw1*t1)*w2*PHI2(j, q)
                                            + h_1s2_east(j)*w1*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q)));
            }
            b(k) = mu.east(j);
            break;
        }
        case BC::Kutta:
        {
            double  w1 = phi1->weightFunction(1);
            double dw1 = phi1->weightFunctionDerivative(1);
            double  w2 = phi2->weightFunction(xi_2(j));
            double dw2 = phi2->weightFunctionDerivative(xi_2(j));
            auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(1, xi_2(j), chi);
            double detJ = dxdx1*dydx2 - dxdx2*dydx1;
            double J11_inv = dydx2/detJ;
            double J21_inv =-dydx1/detJ;
            for (size_t p = 0; p < nx; p++)
            {
                double  t1 = phi1->right(p);
                double dt1 = phi1->rightDerivative(p);
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = J11_inv*(w1*dt1 + dw1*t1)*PHI2(j, q)
                                 + J21_inv*t1*(w2*dPHI2(j, q) + dw2*PHI2(j, q));
            }
            b(k) = mu.east(j);
            break;
        }
        case BC::None:
            break;
    }
}

template void Wing::boundary<int>(const Direction, const BC, const int);
template void Wing::boundary<size_t>(const Direction, const BC, const size_t);
template void Wing::boundary<double>(const Direction, const BC, const double);
template void Wing::boundary<arma::vec>(const Direction, const BC, const arma::vec);
template void Wing::boundary<int>(const Direction, const BC, const double, const double, const int);
template void Wing::boundary<size_t>(const Direction, const BC, const double, const double, const size_t);
template void Wing::boundary<double>(const Direction, const BC, const double, const double, const double);
template void Wing::boundary<arma::vec>(const Direction, const BC, const double, const double, const arma::vec);
template void Wing::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const int);
template void Wing::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const size_t);
template void Wing::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double);
template void Wing::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const arma::vec);
template void Wing::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const int);
template void Wing::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const size_t);
template void Wing::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const double);
template void Wing::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const arma::vec);