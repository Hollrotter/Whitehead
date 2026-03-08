#pragma once
#include "Chebyshev.hpp"

class Airfoil
{
    double c; // Chord length
    double qdyn = 1; // Dynamic pressure
    arma::vec alpha; // Pitch
    size_t nx = 1;
    arma::vec xi = cos(arma::datum::pi*(2*arma::regspace(1, nx) - 1)/(2*nx));
    arma::vec x = c/2*(1 + xi); // x-Coordinates of nodes
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
public:
    Airfoil() = default;
    Airfoil(double _c, size_t _nx) : c(_c), nx(_nx) {}
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