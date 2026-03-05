#pragma once
#include "Camber.hpp"

class VLM
{
    arma::mat x = arma::zeros(2, 2); // x-coordinates of nodes
    arma::vec y = arma::zeros(2);    // y-coordinates of nodes
    size_t nx = x.n_rows-1; // Number of Panels in x-Direction
    size_t ny = y.size()-1; // Number of Panels in y-Direction
    size_t nxy = nx*ny; // Product of nx and ny
    size_t con = 1; // Number of configurations
    arma::mat z = arma::zeros(nx+1, ny+1); // z-coordinates of nodes
    arma::vec z0; // z-coordinates in y-direction (without camber)
    double qdyn = 1; // Dynamic pressure
    arma::vec alpha = arma::zeros(con); // Pitch
    Camber c; // Camber of the surface
    arma::mat A = arma::zeros(nxy, nxy); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::mat nC; // Unit normal vector of the wing
    arma::mat RA = arma::zeros(nxy, 3); // First corner of vortex
    arma::mat RB = arma::zeros(nxy, 3); // Second corner of vortex
    arma::mat RC = arma::zeros(nxy, 3); // Collocation points as vector
    arma::cube rC = arma::zeros(nx, ny, 3); // Collocation points as matrix
    arma::cube rG = arma::zeros(nx, ny, 3); // Pressure point of panel
    arma::cube dcp = arma::zeros(nx+1, ny+1, con); // Difference of non-dimensional pressure
    double area = 0; // Wing area
    arma::vec lift   = arma::zeros(con); // Lift in N
    arma::vec moment = arma::zeros(con); // Moment in Nm
    arma::vec ar = arma::zeros(ny+1, 1); // Rigging (implemented for linear analysis only)
    Symmetry sym = Symmetry::none; // Symmetry (no symmetry or symmetry in y direction)
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
public:
    VLM() = default;
    /// Constructor for an aerodynamic model without z-coordinate
    VLM(arma::mat _x, arma::vec _y) : VLM(_x, _y, arma::zeros(_x.n_cols)) {}
    /// Constructor for an aerodynamic model with z-coordinate
    VLM(arma::mat _x, arma::vec _y, arma::vec _z) : x(_x), y(_y), z0(_z) {}
    // Sets the dynamic pressure
    void dynamicPressure(double);
    // Sets the pitch in degree
    void pitch(double _alpha)
    {
        pitch(_alpha*arma::ones(1));
    }
    // Sets the pitch for multiple configurations
    void pitch(arma::vec);
    // Sets the camber of the wing
    void camber(Camber _c)
    {
        c = _c;
    }
    // Sets the rigging in degree
    void rigging(std::array<double, 2> a)
    {
        ar = arma::datum::pi/180*(a[0] + (a[1]-a[0])*(y-y(0))/(y(ny)-y(0)));
    }
    // Define a symmetry
    void symmetry(Symmetry _sym)
    {
        sym = _sym;
    }
    // Set analysis type to linear or nonlinear
    void aerodynamics(Analysis _analysis)
    {
        analysis = _analysis;
    }
    // Vortex Lattice Method
    void vlm();
    // Gets the lift
    arma::vec get_lift() const
    {
        return lift;
    }
    // Gets the moment
    arma::vec get_moment() const
    {
        return moment;
    }
    // Gets the area
    double get_area() const
    {
        return area;
    }
    // Gets the difference of nondimensional pressure
    arma::cube get_dcp() const
    {
        return dcp;
    }
    // Output y, x and dcp for surface plots
    void output(std::string);
    void operator()(Camber _c)
    {
        c = _c;
    }
    void operator()(Symmetry _sym)
    {
        sym = _sym;
    }
    void operator()(Analysis _analysis)
    {
        analysis = _analysis;
    }
private:
    // Computing some geometrc properties like the collocation points
    void geometry();
    // Calculating the velocity induced by a finite line vortex
    arma::mat line(arma::rowvec, arma::rowvec);
    // Calculating the velocity induced by a horseshoe vortex (used for linear analysis)
    arma::mat horseshoe(arma::rowvec, arma::rowvec);
    // Calculating the velocity induced by a ring vortex (used for nonlinear analysis)
    arma::mat ring(arma::field<arma::rowvec>);
    // Calculating the velocity induced by a semi-infiinte line vortex
    arma::mat wake(size_t, arma::rowvec, arma::rowvec);
    // Calculating the Aerodynamic Matrix for the Vortex-Lattice-Method
    void aerodynamicMatrix();
    // Compute LU-decomposition of system matrix for linear analysis
    void vlmSolve();
    // Compute aerodynamic properties from pre-calculated LU-decomposition
    void vlmEval();
    // Nonlinear Vortex Lattice Method
    void vlmNonlinear();
    // Calculate cL, cM and dcp from gamma
    void postprocessing(arma::mat&);
};