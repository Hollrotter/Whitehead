#include "Aerodynamics.hpp"

template <class C> void Aerodynamics::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const C val)
{
    for (auto &wing:wings)
        for (const auto &chi:wing->chi)
            if (chi == dir)
                wing->boundary(chi, bc, val);
}

template <class C> void Aerodynamics::boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2, const C val)
{
    for (auto &wing:wings)
        for (const auto &chi:wing->chi)
            if (chi == dir)
                wing->boundary(chi, bc, _r1, _r2, val);
}

template void Aerodynamics::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const int);
template void Aerodynamics::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const size_t);
template void Aerodynamics::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double);
template void Aerodynamics::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const arma::vec);
template void Aerodynamics::boundary<int>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const int);
template void Aerodynamics::boundary<size_t>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const size_t);
template void Aerodynamics::boundary<double>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const double);
template void Aerodynamics::boundary<arma::vec>(const Lagrange::CurveInterpolant*, const BC, const double, const double, const arma::vec);