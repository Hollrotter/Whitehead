#include "Structure.hpp"

void Structure::checkMesh()
{
    for (size_t i = 0; i < membranes.size(); i++)
    {
        printf("Checking mesh of membrane number %lu\n", i);
        membranes[i]->checkMesh();
    }
}

template <class C> void Structure::boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const C val)
{
    for (auto &membrane:membranes)
        for (auto &chi:membrane->chi)
            if (chi == dir)
                membrane->boundary(field, chi, bc, val);
}

template <class C> void Structure::boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2, const C val)
{
    for (auto &membrane:membranes)
        for (auto &chi:membrane->chi)
            if (chi == dir)
                membrane->boundary(field, chi, bc, _r1, _r2, val);
}

void Structure::principalStresses(const std::string nxy, const std::string n12)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->principalStresses(nxy+"_"+std::to_string(k), n12+"_"+std::to_string(k));
}

void Structure::principalStrains(const std::string vxy, const std::string gxy, const std::string g12)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->principalStrains(vxy+"_"+std::to_string(k), gxy+"_"+std::to_string(k), g12+"_"+std::to_string(k));
}

void Structure::output(const std::string filename, const Field field)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->output(filename+"_"+std::to_string(k), field);
}

template void Structure::boundary<int>(const Field, const Lagrange::CurveInterpolant*, const BC, const int);
template void Structure::boundary<size_t>(const Field, const Lagrange::CurveInterpolant*, const BC, const size_t);
template void Structure::boundary<double>(const Field, const Lagrange::CurveInterpolant*, const BC, const double);
template void Structure::boundary<arma::vec>(const Field, const Lagrange::CurveInterpolant*, const BC, const arma::vec);
template void Structure::boundary<int>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const int);
template void Structure::boundary<size_t>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const size_t);
template void Structure::boundary<double>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const double);
template void Structure::boundary<arma::vec>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const arma::vec);