#pragma once
#include "Membrane.hpp"
#include "Interface.hpp"

class Structure
{
    std::vector<Membrane*> membranes;
    std::vector<Interface> interfaces;
    size_t iterations = 1000;
    double residualTarget = 1e-10;
    size_t substeps = 1;
    Analysis analysis = Analysis::linear;
    double gamma0 = 2;
public:
    Structure() = default;
    Structure(std::vector<Membrane*> _membranes) : membranes(_membranes)
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
    };
    // Sets the number of substeps for nonlinear analysis (default: substeps = 1)
    void substepControl(const double _substeps)
    {
        substeps = _substeps;
    }
    void setgamma0(double g)
    {
        gamma0 = g;
    }
    void youngsModulus(const double _Et)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->youngsModulus(_Et); } );
    }
    void poissonsRatio(const double _nu)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->poissonsRatio(_nu); } );
    }
    void load(const double f)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->p.fill(f); } );
    }
    void inPlane1(const double _p1)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->p1.fill(_p1); } );
    }
    void inPlane2(const double _p2)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->p1.fill(_p2); } );
    }
    void inPlane11(const double _n11)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->n11.fill(_n11); } );
    }
    void inPlane12(const double _n12)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->n12.fill(_n12); } );
    }
    void inPlane22(const double _n22)
    {
        std::for_each(membranes.begin(), membranes.end(), [&](auto& m) { m->n22.fill(_n22); } );
    }
    void checkMesh();
    void setIterations(const size_t itt)
    {
        iterations = itt;
    }
    void boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc)
    {
        boundary(field, dir, bc, 0);
    }
    template <class C> void boundary(const Field, const Lagrange::CurveInterpolant*, const BC, const C);
    void boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const C);

    void boundary(const Field field, const BC bc, const Lagrange::CurveInterpolant* dir)
    {
        boundary(field, dir, bc, 0);
    }
    template <class C> void boundary(const Field field, const BC bc, const Lagrange::CurveInterpolant* dir, const C c)
    {
        boundary(field, dir, bc, c);
    }
    void boundary(const Field field, const BC bc, const Lagrange::CurveInterpolant* dir, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Field field, const BC bc, const Lagrange::CurveInterpolant* dir, const double _r1, const double _r2, const C c)
    {
        boundary(field, dir, bc, _r1, _r2, c);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const Field field, const BC bc)
    {
        boundary(field, dir, bc, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant* dir, const Field field, const BC bc, const C c)
    {
        boundary(field, dir, bc, c);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const Field field, const BC bc, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant* dir, const Field field, const BC bc, const double _r1, const double _r2, const C c)
    {
        boundary(field, dir, bc, _r1, _r2, c);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const Field field)
    {
        boundary(field, dir, bc, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const Field field, const C c)
    {
        boundary(field, dir, bc, c);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const Field field, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const Field field, const double _r1, const double _r2, const C c)
    {
        boundary(field, dir, bc, _r1, _r2, c);
    }
    void planeStrain();
    void linear();
    void semilinear();
    void nonlinear();
    void principalStresses(const std::string, const std::string);
    void principalStrains(const std::string, const std::string, const std::string);
    void output(const std::string, const Field);
    void output(const Field field, const std::string filename)
    {
        output(filename, field);
    }
};