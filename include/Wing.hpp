#pragma once
#include "TensorField.hpp"
#include "Metric.hpp"
#include "fastgl.h"

class Wing
{
    std::array<Lagrange::CurveInterpolant*, 4> chi;
    arma::mat x = arma::zeros(2, 2); // x-coordinates of nodes
    arma::mat y = arma::zeros(2, 2); // y-coordinates of nodes
    size_t nx = x.n_rows; // Number of nodes in x-Direction
    size_t ny = y.n_cols; // Number of nodes in y-Direction
    size_t nxy = nx*ny; // Product of nx and ny
    TensorField mu{nx, ny}; // Doublet distribution
    arma::vec x1   = Chebyshev::gaussLobatto(nx);
    arma::vec x2   = Chebyshev::gaussLobatto(ny);
    arma::vec xi_1 = Chebyshev::gauss(nx); // Collocation points 1-coordinates
    arma::vec xi_2 = Chebyshev::gauss(ny); // Collocation points 2-coordinates
    arma::mat D1 = Lagrange::derivativeMatrix(xi_1);
    arma::mat D2 = Lagrange::derivativeMatrix(xi_2);
    size_t con = 1; // Number of configurations
    double qdyn = 1; // Dynamic pressure
    arma::vec alpha = arma::zeros(con); // Pitch
    arma::field<arma::mat> J = Jacobian(chi);
    arma::cube e_c = MetricCo(J); // Covariant metric tensor of the surface
    arma::cube ec  = MetricContra(e_c); // Contravariant metric tensor of the surface
    arma::rowvec h_1s1_east = sqrt(ec.slice(0).row(nx-1));
    arma::rowvec h_1s2_east = ec.slice(1).row(nx-1)/sqrt(ec.slice(0).row(nx-1));
    arma::rowvec h_1s1_west = sqrt(ec.slice(0).row(0));
    arma::rowvec h_1s2_west = ec.slice(1).row(0)/sqrt(ec.slice(0).row(0));
    arma::vec h_2s2_south = sqrt(ec.slice(2).col(0));
    arma::vec h_2s1_south = ec.slice(1).col(0)/sqrt(ec.slice(2).col(0));
    arma::vec h_2s2_north = sqrt(ec.slice(2).col(ny-1));
    arma::vec h_2s1_north = ec.slice(1).col(ny-1)/sqrt(ec.slice(2).col(ny-1));
    arma::mat A = arma::zeros(nxy, nxy); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::mat b = arma::zeros(nxy, con);
    arma::mat nC; // Unit normal vector of the wing
    arma::mat mu_hat = arma::zeros(nxy, con); // Amplitudes of doublet distribution
    arma::cube dcp = arma::zeros(nx, ny, con); // Difference of non-dimensional pressure
    double area = 0; // Wing area
    arma::vec lift   = arma::zeros(con); // Lift in N
    arma::vec moment = arma::zeros(con); // Moment in Nm
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
    Wing fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4> _chi)
    {
        auto [_x, _y] = Lagrange::TransfiniteQuadMap(_chi);
        return {_x, _y, _chi};
    }
public:
    Wing() = default;
    Wing(arma::mat _x, arma::mat _y) : x(_x), y(_y) {}
    Wing(arma::mat _x, arma::mat _y, std::array<Lagrange::CurveInterpolant*, 4> _chi) : x(_x), y(_y), chi(_chi) {}
    Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi) : Wing(fromTransfiniteQuadMap(_chi)) {}
    // Sets the dynamic pressure
    void dynamicPressure(double);
    // Sets the pitch in degree
    void pitch(double _alpha)
    {
        pitch(_alpha*arma::ones(1));
    }
    // Sets the pitch for multiple configurations
    void pitch(arma::vec);
    void linear();
    void nonlinear();
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
    // Output x, y and dcp for surface plots
    void output(std::string);
private:
    arma::vec externalContour(double, double, double, double, double, double, arma::vec);
    // Calculating the Aerodynamic Matrix for the Panel Method
    void aerodynamicMatrix();
    void muBoundarySouth(const size_t);
    void muBoundaryNorth(const size_t);
    void muBoundaryWest(const size_t);
    void muBoundaryEast(const size_t);
    void linearSolve();
    void linearEval();
    // Calculate cL, cM and dcp from mu_hat
    void postprocessing();
};