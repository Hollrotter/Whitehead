#pragma once
#include "ChebyshevT.hpp"
#include "ChebyshevU.hpp"
#include "JacobiAlpha.hpp"
#include "JacobiBeta.hpp"
#include "TensorField.hpp"
#include "Metric.hpp"
#include "Wake.hpp"

class Wing
{
    std::array<Lagrange::CurveInterpolant*, 4> chi;
    arma::mat x; // x-coordinates of nodes
    arma::mat y; // y-coordinates of nodes
    size_t nx = x.n_rows; // Number of nodes in x-Direction
    size_t ny = y.n_cols; // Number of nodes in y-Direction
    size_t mx = 5;
    size_t my = 5;
    arma::mat z = arma::zeros(nx, ny); // z-coordinates of nodes
    size_t nxy = nx*ny; // Product of nx and ny
    TensorField mu{nx, ny}; // Doublet distribution
    arma::vec x1 = Chebyshev::gaussLobatto(nx);
    arma::vec x2 = Chebyshev::gaussLobatto(ny);
    std::unique_ptr<BasisFunction> phi1;
    std::unique_ptr<BasisFunction> phi2;
    arma::vec xi_1; // Collocation points 1-coordinates
    arma::vec xi_2; // Collocation points 2-coordinates
    arma::mat  PHI1 = arma::mat(nx, nx, arma::fill::none);
    arma::mat  PHI2 = arma::mat(ny, ny, arma::fill::none);
    arma::mat dPHI1 = arma::mat(nx, nx, arma::fill::none);
    arma::mat dPHI2 = arma::mat(ny, ny, arma::fill::none);
    arma::mat D1;
    arma::mat D2;
    double qdyn = 1; // Dynamic pressure
    double alpha = 0; // Pitch
    static const size_t n_theta = 97;
    static constexpr double delta = 0.5;
    arma::mat xC;
    arma::mat yC;
    arma::mat Tx;
    arma::mat Ty;
    arma::mat zC;
    arma::field<arma::mat> J = Jacobian(x1, x2, chi);
    arma::cube e_c = MetricCo(J); // Covariant metric tensor of the surface
    arma::cube ec  = MetricContra(e_c); // Contravariant metric tensor of the surface
    arma::vec    h_2s2_south;
    arma::vec    h_2s1_south;
    arma::rowvec h_1s1_east;
    arma::rowvec h_1s2_east;
    arma::vec    h_2s2_north;
    arma::vec    h_2s1_north;
    arma::rowvec h_1s1_west;
    arma::rowvec h_1s2_west;
    arma::mat A = arma::zeros(nxy, nxy); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::vec b = arma::zeros(nxy);
    arma::mat nC; // Unit normal vector of the wing
    arma::vec mu_hat = arma::zeros(nxy); // Amplitudes of doublet distribution
    arma::mat dcp = arma::zeros(nx, ny); // Difference of non-dimensional pressure
    double area   = 0; // Wing area
    double lift   = 0; // Lift in N
    double moment = 0; // Moment in Nm
    Symmetry sym = Symmetry::none; // Symmetry (no symmetry or symmetry in the y-direction)
    Analysis analysis = Analysis::linear; // Analysis type (linear or nonlinear)
    std::vector<Wake*> wakes;
    Wing fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4>);
    Wing fromTransfiniteQuadMap(arma::mat, std::array<Lagrange::CurveInterpolant*, 4>);
public:
    Wing(arma::mat _x, arma::mat _y) : x(_x), y(_y) {}
    Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi, arma::mat _x, arma::mat _y) : chi(_chi), x(_x), y(_y) {}
    Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi, arma::mat _x, arma::mat _y, arma::mat _z) : chi(_chi), x(_x), y(_y), z(_z) {}
    explicit Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi) : Wing(fromTransfiniteQuadMap(_chi)) {}
    Wing(arma::mat _z, std::array<Lagrange::CurveInterpolant*, 4> _chi) : Wing(fromTransfiniteQuadMap(_z, _chi)) {}
    // Sets the dynamic pressure
    void dynamicPressure(double _qdyn) pre(_qdyn > 0 && "Dynamic pressure must be positive!")
    {
        qdyn = _qdyn;
    }
    // Sets the pitch in degree
    void pitch(double _alpha)
    {
        alpha = arma::datum::pi/180*_alpha;
    }
    // Define a symmetry
    void symmetry(Symmetry _sym)
    {
        sym = _sym;
    }
    void checkMesh() const;
    void wake(Wake* w)
    {
        wakes.push_back(w);
    }
    void linear();
    void nonlinear();
    // Gets the lift
    double get_lift() const
    {
        return lift;
    }
    // Gets the moment
    double get_moment() const
    {
        return moment;
    }
    // Gets the area
    double get_area() const
    {
        return area;
    }
    // Gets the difference of nondimensional pressure
    arma::mat get_dcp() const
    {
        return dcp;
    }
    void boundary(const Direction dir, const BC bc)
    {
        boundary(dir, bc, 0);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc)
    {
        boundary(dir, bc, 0);
    }
    template <class C> void boundary(const Direction, const BC, const C);
    template <class C> void boundary(const Lagrange::CurveInterpolant*, const BC, const C);
    void boundary(const Direction dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(dir, bc, _r1, _r2, 0);
    }
    void boundary(const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(dir, bc, _r1, _r2, 0);
    }
    template <class C> void boundary(const Direction, const BC, const double, const double, const C);
    template <class C> void boundary(const Lagrange::CurveInterpolant*, const BC, const double, const double, const C);
    // Output x, y and dcp for surface plots
    void output(std::string) const;
    void operator()(Symmetry _sym)
    {
        sym = _sym;
    }
    friend class Aerodynamics;
private:
    void init();
    arma::mat calculateNormal();
    arma::vec externalContour(double, double, double, double, double, double, arma::vec);
    void regularIntegralLinear(size_t, double, double, size_t, size_t, double, double, double, double);
    void regularIntegralNonlinear(size_t, double, double, double, size_t, size_t, double, double, double, double);
    // Calculating the Aerodynamic Matrix for the Panel Method
    void aerodynamicMatrix();
    void muBoundarySouth(const size_t);
    void muBoundaryNorth(const size_t);
    void muBoundaryWest(const size_t);
    void muBoundaryEast(const size_t);
    void linearSolve();
    void linearEval();
    void nonlinearSolve();
    void nonlinearEval();
    // Calculate cL, cM and dcp from mu_hat
    void postprocessing();
};