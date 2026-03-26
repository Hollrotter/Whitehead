#pragma once
#include "Lagrange.hpp"

class Airfoil
{
    double c; // Chord length
    double qdyn = 1; // Dynamic pressure
    arma::vec alpha; // Pitch
    size_t nx = 1;
    Lagrange::CurveInterpolant* chi;
    arma::vec xi = Chebyshev::gauss(nx);
    arma::vec x = c/2*(1 + xi); // x-Coordinates of nodes
    arma::vec z = arma::zeros(nx);
    size_t con = 1; // Number of configurations
    arma::mat gamma_hat;
    arma::mat dcp; // Difference of non-dimensional pressure
    arma::vec cL; // Lift coefficient
    arma::vec cM; // Moment coefficient
    arma::mat A = arma::zeros(nx, nx); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::mat b = arma::zeros(nx, con);
    arma::mat nC; // Normal vector of the airfoil
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
    Airfoil fromLagrangeCurveInterpolant(Lagrange::CurveInterpolant* _chi)
    {
        arma::vec _xi = Chebyshev::gauss(_chi->getNodes().size());
        auto [_x, _z] = _chi->evaluate(_xi);
        return {_x, _z, _chi};
    }
public:
    Airfoil() = default;
    Airfoil(double _c, size_t _nx) : c(_c), nx(_nx) {}
    Airfoil(arma::vec _x, arma::vec _z, Lagrange::CurveInterpolant* _chi) : x(_x), z(_z), chi(_chi), c(_x.back()-_x.front()), nx(_x.size()) {};
    Airfoil(Lagrange::CurveInterpolant* _chi) : Airfoil(fromLagrangeCurveInterpolant(_chi)) {}
    // Set dynamic pressure
    void dynamicPressure(double);
    // Set pitch in degree
    void pitch(double);
    // Set pitch in degree
    void pitch(arma::vec);
    void linear();
    void nonlinear();
    arma::vec get_lift() const
    {
        return cL;
    }
    arma::vec get_moment() const
    {
        return cM;
    }
    arma::vec get_dcp() const
    {
        return dcp;
    }
    // output x and dcp to given file
    void output(std::string);
private:
    // Calculates the Aerodynamic Matrix needed for the Discrete-Vortex-Method
    void aerodynamicMatrix();
    void linearSolve();
    void linearEval();
    void postprocessing();
};