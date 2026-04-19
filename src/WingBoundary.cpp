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
            for (size_t p = 0; p < nx; p++)
            {
                double T = boost::math::chebyshev_t(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(i, p+q*nx) = T * pow(-1, q);
            }
            b(i) = mu.south(i);
            break;
        case BC::Neumann:
            for (size_t p = 0; p < nx; p++)
            {
                double  T = boost::math::chebyshev_t(p, xi_1(i));
                double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(i, p+q*nx) = h_2s1_south(i)*dT*pow(-1, q) + h_2s2_south(i)*T*pow(-1, q+1)*pow(q, 2);
            }
            b(i) = mu.south(i);
            break;
        case BC::Robin:
            for (size_t p = 0; p < nx; p++)
            {
                double  T = boost::math::chebyshev_t(p, xi_1(i));
                double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(i, p+q*nx) = mu.r1South*T*pow(-1, q)
                                 + mu.r2South*(h_2s1_south(i)*dT*pow(-1, q) + h_2s2_south(i)*T*pow(-1, q+1)*pow(q, 2));
            }
            b(i) = mu.south(i);
            break;
        case BC::Derivative_x:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), -1, chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J11_inv = dydx2/detJ;
                double J21_inv =-dydx1/detJ;
                for (size_t p = 0; p < nx; p++)
                {
                    double  T = boost::math::chebyshev_t(p, xi_1(i));
                    double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                    for (size_t q = 0; q < ny; q++)
                        A(i, p+q*nx) = J11_inv*dT*pow(-1, q) + J21_inv*T*pow(-1, q+1)*pow(q, 2);
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(i) = mu.south(i);
            break;
        case BC::Derivative_y:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), -1, chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J12_inv =-dxdx2/detJ;
                double J22_inv = dxdx1/detJ;
                for (size_t p = 0; p < nx; p++)
                {
                    double  T = boost::math::chebyshev_t(p, xi_1(i));
                    double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                    for (size_t q = 0; q < ny; q++)
                        A(i, p+q*nx) = J12_inv*dT*pow(-1, q) + J22_inv*T*pow(-1, q+1)*pow(q, 2);
                }
            }
            else
            {
                std::println("Derivative_y not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(i) = mu.south(i);
            break;
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
            for (size_t p = 0; p < nx; p++)
            {
                double T = boost::math::chebyshev_t(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = T;
            }
            b(k) = mu.north(i);
            break;
        case BC::Neumann:
            for (size_t p = 0; p < nx; p++)
            {
                double  T = boost::math::chebyshev_t(p, xi_1(i));
                double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = h_2s1_north(i)*dT + h_2s2_north(i)*T*pow(q, 2);
            }
            b(k) = mu.north(i);
            break;
        case BC::Robin:
            for (size_t p = 0; p < nx; p++)
            {
                double  T = boost::math::chebyshev_t(p, xi_1(i));
                double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                for (size_t q = 0; q < ny; q++)
                    A(k, p+q*nx) = mu.r1North*T + mu.r2North*(h_2s1_north(i)*dT + h_2s2_north(i)*T*pow(q, 2));
            }
            b(k) = mu.north(i);
            break;
        case BC::Derivative_x:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), 1, chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J11_inv = dydx2/detJ;
                double J21_inv =-dydx1/detJ;
                for (size_t p = 0; p < nx; p++)
                {
                    double  T = boost::math::chebyshev_t(p, xi_1(i));
                    double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                    for (size_t q = 0; q < ny; q++)
                        A(k, p+q*nx) = J11_inv*dT + J21_inv*T*pow(q, 2);
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.north(i);
            break;
        case BC::Derivative_y:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(xi_1(i), 1, chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J12_inv =-dxdx2/detJ;
                double J22_inv = dxdx1/detJ;
                for (size_t p = 0; p < nx; p++)
                {
                    double  T = boost::math::chebyshev_t(p, xi_1(i));
                    double dT = boost::math::chebyshev_t_prime(p, xi_1(i));
                    for (size_t q = 0; q < ny; q++)
                        A(k, p+q*nx) = J12_inv*dT + J22_inv*T*pow(q, 2);
                }
            }
            else
            {
                std::println("Derivative_y not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.north(i);
            break;
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
            for (size_t q = 0; q < ny; q++)
            {
                double T = boost::math::chebyshev_t(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = pow(-1, p) * T;
            }
            b(k) = mu.west(j);
            break;
        case BC::Neumann:
            for (size_t q = 0; q < ny; q++)
            {
                double  T = boost::math::chebyshev_t(q, xi_2(j));
                double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = h_1s1_west(j)*pow(-1, p+1)*pow(p, 2)*T + h_1s2_west(j)*pow(-1, p)*dT;
            }
            b(k) = mu.west(j);
            break;
        case BC::Robin:
            for (size_t q = 0; q < ny; q++)
            {
                double  T = boost::math::chebyshev_t(q, xi_2(j));
                double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = mu.r1West*pow(-1, p) * T
                                 + mu.r2West*(h_1s1_west(j)*pow(-1, p+1)*pow(p, 2)*T + h_1s2_west(j)*pow(-1, p)*dT);
            }
            b(k) = mu.west(j);
            break;
        case BC::Derivative_x:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(-1, xi_2(j), chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J11_inv = dydx2/detJ;
                double J21_inv =-dydx1/detJ;
                for (size_t q = 0; q < ny; q++)
                {
                    double  T = boost::math::chebyshev_t(q, xi_2(j));
                    double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                    for (size_t p = 0; p < nx; p++)
                        A(k, p+q*nx) = J11_inv*pow(-1, p+1)*pow(p, 2)*T + J21_inv*pow(-1, p)*dT;
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.west(j);
            break;
        case BC::Derivative_y:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(-1, xi_2(j), chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J12_inv =-dxdx2/detJ;
                double J22_inv = dxdx1/detJ;
                for (size_t q = 0; q < ny; q++)
                {
                    double  T = boost::math::chebyshev_t(q, xi_2(j));
                    double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                    for (size_t p = 0; p < nx; p++)
                        A(k, p+q*nx) = J12_inv*pow(-1, p+1)*pow(p, 2)*T + J22_inv*pow(-1, p)*dT;
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.west(j);
            break;
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
            for (size_t q = 0; q < ny; q++)
            {
                double T = boost::math::chebyshev_t(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = T;
            }
            b(k) = mu.east(j);
            break;
        case BC::Neumann:
            for (size_t q = 0; q < ny; q++)
            {
                double  T = boost::math::chebyshev_t(q, xi_2(j));
                double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = h_1s1_east(j)*pow(p, 2)*T + h_1s2_east(j)*dT;
            }
            b(k) = mu.east(j);
            break;
        case BC::Robin:
            for (size_t q = 0; q < ny; q++)
            {
                double  T = boost::math::chebyshev_t(q, xi_2(j));
                double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                for (size_t p = 0; p < nx; p++)
                    A(k, p+q*nx) = mu.r1East*T + mu.r2East*(h_1s1_east(j)*pow(p, 2)*T + h_1s2_east(j)*dT);
            }
            b(k) = mu.east(j);
            break;
        case BC::Derivative_x:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(1, xi_2(j), chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J11_inv = dydx2/detJ;
                double J21_inv =-dydx1/detJ;
                for (size_t q = 0; q < ny; q++)
                {
                    double  T = boost::math::chebyshev_t(q, xi_2(j));
                    double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                    for (size_t p = 0; p < nx; p++)
                        A(k, p+q*nx) = J11_inv*pow(p, 2)*T + J21_inv*dT;
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.east(j);
            break;
        case BC::Derivative_y:
            if (analysis == Analysis::linear)
            {
                auto [dxdx1, dxdx2, dydx1, dydx2] = Lagrange::TransfiniteQuadMetrics(1, xi_2(j), chi);
                double detJ = dxdx1*dydx2 - dxdx2*dydx1;
                double J12_inv =-dxdx2/detJ;
                double J22_inv = dxdx1/detJ;
                for (size_t q = 0; q < ny; q++)
                {
                    double  T = boost::math::chebyshev_t(q, xi_2(j));
                    double dT = boost::math::chebyshev_t_prime(q, xi_2(j));
                    for (size_t p = 0; p < nx; p++)
                        A(k, p+q*nx) = J12_inv*pow(p, 2)*T + J22_inv*dT;
                }
            }
            else
            {
                std::println("Derivative_x not yet implemented for nonlinear analysis!");
                exit(EXIT_FAILURE);
            }
            b(k) = mu.east(j);
            break;
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