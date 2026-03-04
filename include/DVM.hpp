#pragma once
#include "Splinefit.hpp"
#include "Chebyshev.hpp"

class DVM
{
    CamberType camberType; //Camber type (none, function or B_Spline)
    std::function<double(double)>  f; // Analytical function for the camber
    std::function<double(double)> df; // Analytical function for the camber slope
    Splinefit s; // Splinefitting for the initial contour of the camber
    double c; // Chord length
    double qdyn = 1; // Dynamic pressure
    arma::vec alpha = arma::zeros(1); // Pitch
    size_t nx = 1; // Number of segments
    size_t con = 1; // Number of configurations
    arma::vec x = c/2*(1 + Chebyshev::gaussLobatto(nx+1)); // x-Coordinates of nodes
    arma::vec z = arma::zeros(nx+1); // z-Coordinates of nodes
    arma::vec xg; // x-coordinates of Pressure points
    arma::vec xC; // x-coordinates of Collocation points
    arma::vec zg; // z-coordinates of Pressure points
    arma::vec zC; // z-coordinates of Collocation points
    arma::mat dcp; // Difference of non-dimensional pressure
    arma::vec cL; // Lift coefficient
    arma::vec cM; // Moment coefficient
    arma::mat A = arma::zeros(nx, nx); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::mat nC; // Normal vector of the surface
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
public:
    DVM() : camberType(CamberType::none) {}
    // Constructor for flat plate
    DVM(double _c, size_t _nx) : c(_c), nx(_nx), camberType(CamberType::none) {}
    // Constructor for camber defined by function describing only the derivative of z (use for linear analysis only!)
    DVM(std::function<double(double)> dF) : df(dF), camberType(CamberType::function) {}
    // Constructor for camber defined by function describing only the derivative of z (use for linear analysis only!)
    DVM(double _c, size_t _nx, std::function<double(double)> dF) : c(_c), nx(_nx), df(dF), camberType(CamberType::function) {}
    // Constructor for camber defined by function describing the value and the derivative of z (use for nonlinear analysis!)
    DVM(std::function<double(double)> F, std::function<double(double)> dF) : f(F), df(dF), camberType(CamberType::function) {}
    // Constructor for camber defined by function describing the value and the derivative of z (use for nonlinear analysis!)
    DVM(double _c, size_t _nx, std::function<double(double)> F, std::function<double(double)> dF) : c(_c), nx(_nx), f(F), df(dF), camberType(CamberType::function) {}
    // Constructor for camber defined by Splinefitting
    DVM(Splinefit S) : s(S), camberType(CamberType::b_spline) {}
    // Constructor for camber defined by Splinefitting
    DVM(double _c, size_t _nx, Splinefit S) : c(_c), nx(_nx), s(S), camberType(CamberType::b_spline) {}
    // Compute z coordinate at given x value
    double operator() (double);
    // Compute z coordinates at given vector x
    arma::vec operator() (arma::vec);
    // Compute derivative at given x value
    double diff(double);
    // Compute derivatives at given vector x
    arma::vec diff(arma::vec);
    // Set dynamic pressure
    void dynamicPressure(double);
    // Set pitch in degree
    void pitch(double);
    // Set pitch in degree
    void pitch(arma::vec);
    // Set analysis type (linear or nonlinear)
    void aerodynamics(Analysis _analysis)
    {
        analysis = _analysis;
    }
    // Discrete Vortex Method
    void dvm();
    // Getter for the lift coefficient
    arma::vec get_lift() const
    {
        return cL;
    }
    // Getter for the moment coefficient
    arma::vec get_moment() const
    {
        return cM;
    }
    // Getter for difference of non-dimensional pressure
    arma::mat get_dcp() const
    {
        return dcp;
    }
    // output x and dcp to given file
    void output(std::string);
    void operator()(Analysis _analysis)
    {
        analysis = _analysis;
    }
private:
    // Computes some geometric properties like the collocation points
    void geometry();
    // Calculates the Aerodynamic Matrix needed for the Discrete-Vortex-Method
    void aerodynamicMatrix();
    // Calculate LU-decomposition for lienar Discrete Vortex Method
    void dvmSolve();
    // Compute results from pre-calculated LU-decomposition
    void dvmEval();
    // Nonlinear Discrete Vortex Method
    void dvmNonlinear();
    // Calculate cL, cM and dcp from gamma
    void postprocessing(arma::mat&);
};