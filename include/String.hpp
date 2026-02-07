#pragma once
#include "Chebyshev.hpp"

class String
{
    double c = 1; // Chord length
    size_t n = 2; // Number of nodes in x-direction
    double sigma = 0; // Pretension
    double Et = 10; // Young's modulus multiplied by thickness
    arma::vec x   = c/2*(Chebyshev::gaussLobatto(n) + 1); // x-coordinates of the nodes
    arma::mat D1  = Chebyshev::derivativeMatrix(Chebyshev::gaussLobatto(n), Derivative::first)*2/c;
    arma::mat D11 = Chebyshev::derivativeMatrix(Chebyshev::gaussLobatto(n), Derivative::second)*4/c/c;
    arma::vec z = arma::zeros(n); // z-coordinates of the nodes
    double L0 = sqrt(c*c + pow(z(n-1)-z(0), 2)); // Unstrained length of the string
    arma::vec p = arma::zeros(n); // Pressure at the nodes
    arma::mat S = arma::zeros(n, n); // Structural Matrix
    arma::mat L = arma::zeros(n, n); // Lower triangular matrix
    arma::mat U = arma::zeros(n, n); // Upper triangular matrix
    arma::mat P = arma::zeros(n, n); // Permutation matrix
    arma::vec b = arma::zeros(n); // Right side of linear equation system
    BC front = BC::Dirichlet; // Boundary type of front
    BC back  = BC::Dirichlet; // Boundary type of back
    double frontBC = 0; // Boundary value at front
    double  backBC = 0; // Boundary value at back
    size_t iter = 150; // Max. # of iterations for nonlinear solution
    double omega = 1; // Underrelaxation for Newton-Method
    double residualTarget = 1e-5; // Target residual for nonlinear solution
    Material materialModel = Material::inextensible;
    Analysis analysis = Analysis::linear;
public:
    String() = default;
    // Constructor for a String giving the chord length, the number of nodes and pretension
    String(double, size_t, double);
    // Sets the young's modulus times thickness (Et)
    void youngsModulus(double _Et);
    // Sets the material model (inextensible or extensible)
    void material(Material _materialModel)
    {
        materialModel = _materialModel;
    }
    // Sets the under-relaxation factor for Newton's Method (default: omega = 1)
    void relaxationFactor(double _omega);
    // Applying a constant load on the string
    void load(double f)
    {
        p.fill(f);
    }
    // Applying a distributed load on the string
    void load(arma::vec f)
    {
        p = f;
    }
    // Applying a distributed load defined by a lambda function
    void load(std::function<double(double)>);
    // Setting a homogenious boundary condition at given location and boundary type
    void boundary(Location loc, BC bc)
    {
        boundary(loc, bc, 0);
    }
    // Setting a boundary condition at given location and boundary type with scalar value
    void boundary(Location, BC, double);
    // Linear analysis of string computing the deformation z
    void linear();
    // Nonlinear analysis of string computing the deformation z
    void nonlinear();
    // Calculates the integral of a given vector over the length of the string
    double integrate(arma::vec);
    // Output x, y and z to a given file
    void output(std::string);
    void operator()(Material _materialModel)
    {
        materialModel = _materialModel;
    }
private:
    void structuralMatrix();
    void solve_S();
    void solve_b();
};