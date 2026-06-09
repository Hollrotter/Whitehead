#pragma once
#include "Wing.hpp"
#include "Interface.hpp"

class Aerodynamics
{
    std::vector<Wing*> wings;
    std::vector<Interface> interfaces;
    size_t iterations = 100;
    double residualTarget = 1e-10;
    Symmetry sym = Symmetry::none;
    Analysis analysis = Analysis::linear;
    double lambda0 = 2;
    arma::field<arma::cube> bw;
    Aerodynamics fromWings(std::vector<Wing*>);
public:
    Aerodynamics() = default;
    Aerodynamics(const std::vector<Wing*> _w, const std::vector<Interface> _i) : wings(_w), interfaces(_i) {};
    Aerodynamics(const std::vector<Wing*> _w) : Aerodynamics(fromWings(_w)) {};
    void setlambda(double);
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
    // Define a symmetry
    void symmetry(Symmetry _sym)
    {
        sym = _sym;
    }
    void checkMesh();
    void wake(Wake*);
    void setIterations(const size_t itt)
    {
        iterations = itt;
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc)
    {
        boundary(dir, bc, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant*, const BC, const C);
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Lagrange::CurveInterpolant*, const BC, const double, const double, const C);
    void linear();
    void nonlinear();
    arma::vec get_lift();
    arma::vec get_moment();
    double get_area();
    void output(const std::string);
    void operator()(Symmetry _sym)
    {
        sym = _sym;
    }
private:
    void solve();
};