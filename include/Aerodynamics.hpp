#pragma once
#include "Wing.hpp"
#include "Interface.hpp"

class Aerodynamics
{
    std::vector<Wing*> wings;
    std::vector<Interface> interfaces;
    size_t iterations = 1000;
    double residualTarget = 1e-10;
    Analysis analysis = Analysis::linear;
    double gamma0 = 2;
    arma::field<arma::cube> bw;
public:
    Aerodynamics() = default;
    Aerodynamics(std::vector<Wing*> _wings) : wings(_wings)
    {
        for (size_t sD = 0; sD < wings.size()-1; sD++)
            for (size_t sC = 0; sC < 4; sC++)
                for (size_t tD = sD+1; tD < wings.size(); tD++)
                    for (size_t tC = 0; tC < 4 ; tC++)
                        if (wings[sD]->chi[sC] == wings[tD]->chi[tC])
                        {
                            interfaces.push_back(Interface(sD, tD, sC, tC));
                            wings[sD]->chi[sC]->curveType = CurveType::Interface;
                        }
    }
    void setgamma0(double g)
    {
        gamma0 = g;
    }
    // Sets the dynamic pressure
    void dynamicPressure(double _qdyn)
    {
        std::for_each(wings.begin(), wings.end(), [&](auto& w) { w->dynamicPressure(_qdyn); } );
    }
    // Sets the pitch in degree
    void pitch(const double _alpha)
    {
        std::for_each(wings.begin(), wings.end(), [&](auto& w) { w->pitch(_alpha); } );
    }
    // Sets the pitch for multiple configurations
    void pitch(const arma::vec _alpha)
    {
        std::for_each(wings.begin(), wings.end(), [&](auto& w) { w->pitch(_alpha); } );
    }
    void checkMesh();
    void setIterations(const size_t itt)
    {
        iterations = itt;
    }
    void linear();
    void nonlinear();
    void output(const std::string);
private:
    void solve();
};