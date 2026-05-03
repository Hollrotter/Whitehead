#include "Structure.hpp"

Structure::Structure(std::vector<Membrane*> _membranes) : membranes(_membranes)
{
    for (size_t sD = 0; sD < membranes.size()-1; sD++)
        for (size_t sC = 0; sC < 4; sC++)
            for (size_t tD = sD+1; tD < membranes.size(); tD++)
                for (size_t tC = 0; tC < 4; tC++)
                    if (membranes[sD]->chi[sC] == membranes[tD]->chi[tC])
                    {
                        interfaces.push_back(Interface(sD, tD, sC, tC));
                        membranes[sD]->chi[sC]->curveType = CurveType::Interface;
                    }
    for (Interface& interface:interfaces)
    {
        Membrane *mSource = membranes[interface.sourceDomain];
        switch (interface.sourceCurve)
        {
            case 0: // South
                interface.c1Source  = mSource->h_1s1_south%mSource->h_2s1_south;
                interface.c2Source  = mSource->h_1s1_south%mSource->h_2s2_south;
                interface.c11Source = pow(mSource->h_2s1_south, 2);
                interface.c12Source = mSource->h_2s1_south%mSource->h_2s2_south;
                interface.c22Source = pow(mSource->h_2s2_south, 2);
                interface.lambdaSource = lambda0*mean(mSource->h_2s2_south);
                break;
            case 1: // East
                interface.c1Source  = (mSource->h_2s2_east%mSource->h_1s1_east).t();
                interface.c2Source  = (mSource->h_2s2_east%mSource->h_1s2_east).t();
                interface.c11Source = pow(mSource->h_1s1_east, 2).t();
                interface.c12Source = (mSource->h_1s1_east%mSource->h_1s2_east).t();
                interface.c22Source = pow(mSource->h_1s2_east, 2).t();
                interface.lambdaSource = lambda0*mean(mSource->h_1s1_east);
                break;
            case 2: // North
                interface.c1Source  = mSource->h_1s1_north%mSource->h_2s1_north;
                interface.c2Source  = mSource->h_1s1_north%mSource->h_2s2_north;
                interface.c11Source = pow(mSource->h_2s1_north, 2);
                interface.c12Source = mSource->h_2s1_north%mSource->h_2s2_north;
                interface.c22Source = pow(mSource->h_2s2_north, 2);
                interface.lambdaSource = lambda0*mean(mSource->h_2s2_north);
                break;
            case 3: // West
                interface.c1Source  = (mSource->h_2s2_west%mSource->h_1s1_west).t();
                interface.c2Source  = (mSource->h_2s2_west%mSource->h_1s2_west).t();
                interface.c11Source = pow(mSource->h_1s1_west, 2).t();
                interface.c12Source = (mSource->h_1s1_west%mSource->h_1s2_west).t();
                interface.c22Source = pow(mSource->h_1s2_west, 2).t();
                interface.lambdaSource = lambda0*mean(mSource->h_1s1_west);
                break;
        }
        Membrane *mTarget = membranes[interface.targetDomain];
        switch (interface.targetCurve)
        {
            case 0: // South
                interface.c1Target  = mTarget->h_1s1_south%mTarget->h_2s1_south;
                interface.c2Target  = mTarget->h_1s1_south%mTarget->h_2s2_south;
                interface.c11Target = pow(mTarget->h_2s1_south, 2);
                interface.c12Target = mTarget->h_2s1_south%mTarget->h_2s2_south;
                interface.c22Target = pow(mTarget->h_2s2_south, 2);
                interface.lambdaTarget = lambda0*mean(mTarget->h_2s2_south);
                break;
            case 1: // East
                interface.c1Target  = (mTarget->h_2s2_east%mTarget->h_1s1_east).t();
                interface.c2Target  = (mTarget->h_2s2_east%mTarget->h_1s2_east).t();
                interface.c11Target = pow(mTarget->h_1s1_east, 2).t();
                interface.c12Target = (mTarget->h_1s1_east%mTarget->h_1s2_east).t();
                interface.c22Target = pow(mTarget->h_1s2_east, 2).t();
                interface.lambdaTarget = lambda0*mean(mTarget->h_1s1_east);
                break;
            case 2: // North
                interface.c1Target  = mTarget->h_1s1_north%mTarget->h_2s1_north;
                interface.c2Target  = mTarget->h_1s1_north%mTarget->h_2s2_north;
                interface.c11Target = pow(mTarget->h_2s1_north, 2);
                interface.c12Target = mTarget->h_2s1_north%mTarget->h_2s2_north;
                interface.c22Target = pow(mTarget->h_2s2_north, 2);
                interface.lambdaTarget = lambda0*mean(mTarget->h_2s2_north);
                break;
            case 3: // West
                interface.c1Target  = (mTarget->h_2s2_west%mTarget->h_1s1_west).t();
                interface.c2Target  = (mTarget->h_2s2_west%mTarget->h_1s2_west).t();
                interface.c11Target = pow(mTarget->h_1s1_west, 2).t();
                interface.c12Target = (mTarget->h_1s1_west%mTarget->h_1s2_west).t();
                interface.c22Target = pow(mTarget->h_1s2_west, 2).t();
                interface.lambdaTarget = lambda0*mean(mTarget->h_1s1_west);
                break;
        }
    }
};

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