#pragma once
#include "Camber.hpp"
#include "Lagrange.hpp"

class Airfoil
{
    Camber camber;
    double c; // Chord length
    double qdyn = 1; // Dynamic pressure
    double alpha; // Pitch
    size_t nx = 1;
    Lagrange::CurveInterpolant* chi;
    arma::vec xi = Chebyshev::gauss(nx);
    arma::vec x = c/2*(1 + xi); // x-Coordinates of nodes
    arma::vec z = arma::zeros(nx);
    arma::vec gamma_hat;
    arma::vec dcp; // Difference of non-dimensional pressure
    double cL; // Lift coefficient
    double cM; // Moment coefficient
    arma::mat A = arma::zeros(nx, nx); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::vec b = arma::zeros(nx);
    arma::mat nC; // Normal vector of the airfoil
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
    Airfoil fromLagrangeCurveInterpolant(Lagrange::CurveInterpolant*);
public:
    Airfoil() : camber(Camber()) {}
    Airfoil(double _c, size_t _nx) : camber(Camber()), c(_c), nx(_nx) {}
    Airfoil(double _c, size_t _nx, std::function<double(double)> dF) : camber(Camber(dF)), c(_c), nx(_nx) {}
    Airfoil(arma::vec _x, arma::vec _z, Lagrange::CurveInterpolant* _chi) : c(_x.back()-_x.front()), nx(_x.size()), chi(_chi), x(_x), z(_z) {};
    explicit Airfoil(Lagrange::CurveInterpolant* _chi) : Airfoil(fromLagrangeCurveInterpolant(_chi)) {}
    // Set dynamic pressure
    void dynamicPressure(double _qdyn) pre(_qdyn > 0 && "Dynamic pressure must be positive!")
    {
        qdyn = _qdyn;
    }
    // Set pitch in degree
    void pitch(double _alpha)
    {
        alpha = arma::datum::pi/180*_alpha;
    }
    void linear();
    void nonlinear();
    double get_lift() const
    {
        return cL;
    }
    double get_moment() const
    {
        return cM;
    }
    arma::vec get_dcp() const
    {
        return dcp;
    }
    // output x and dcp to given file
    void output(std::string) const;
private:
    double x1(double, double) const;
    double z1(double, double) const;
    double x2(double, double, double, double) const;
    double z2(double, double, double, double) const;
    double x3(double, double, double, double, double) const;
    double z3(double, double, double, double, double) const;
    double r2(double, double, double, double, double, double) const;
    double k1(double, double, double, double, double, double, double, double, double, double) const;
    double k2(double, double, double, double, double, double, double, double, double, double) const;
    double k3(double, double, double, double, double, double, double, double, double) const;
    // Calculates the Aerodynamic Matrix needed for the Discrete-Vortex-Method
    void aerodynamicMatrix();
    void linearSolve();
    void linearEval();
    void postprocessing();
};