#pragma once
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
    arma::mat z = arma::zeros(nx, ny); // z-coordinates of nodes
    size_t nxy = nx*ny; // Product of nx and ny
    TensorField mu{nx, ny}; // Doublet distribution
    arma::vec x1   = Chebyshev::gaussLobatto(nx);
    arma::vec x2   = Chebyshev::gaussLobatto(ny);
    arma::vec xi_1 = Chebyshev::gauss(nx); // Collocation points 1-coordinates
    arma::vec xi_2 = Chebyshev::gauss(ny); // Collocation points 2-coordinates
    arma::mat  T1 = Chebyshev::Polynomial(xi_1);
    arma::mat  T2 = Chebyshev::Polynomial(xi_2);
    arma::mat dT1 = Chebyshev::derivative(xi_1);
    arma::mat dT2 = Chebyshev::derivative(xi_2);
    arma::mat D1 = Lagrange::derivativeMatrix(xi_1);
    arma::mat D2 = Lagrange::derivativeMatrix(xi_2);
    double qdyn = 1; // Dynamic pressure
    double alpha = 0; // Pitch
    size_t n_theta = 50;
    size_t m_rho = 2;
    std::tuple<arma::mat, arma::mat> xyC = Lagrange::TransfiniteQuadMap(xi_1, xi_2, chi);
    arma::mat xC = std::get<0>(xyC);
    arma::mat yC = std::get<1>(xyC);
    arma::mat Tx = Lagrange::interpolationMatrix(x1, xi_1);
    arma::mat Ty = Lagrange::interpolationMatrix(x2, xi_2);
    arma::mat zC = Lagrange::interpolation2D(Tx, Ty, z, xi_1, xi_2);
    arma::field<arma::mat> J = Jacobian(x1, x2, chi);
    arma::cube e_c = MetricCo(J); // Covariant metric tensor of the surface
    arma::cube ec  = MetricContra(e_c); // Contravariant metric tensor of the surface
    std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec> h;
    arma::vec    h_2s2_south = std::get<0>(h);
    arma::vec    h_2s1_south = std::get<1>(h);
    arma::rowvec h_1s1_east  = std::get<2>(h);
    arma::rowvec h_1s2_east  = std::get<3>(h);
    arma::vec    h_2s2_north = std::get<4>(h);
    arma::vec    h_2s1_north = std::get<5>(h);
    arma::rowvec h_1s1_west  = std::get<6>(h);
    arma::rowvec h_1s2_west  = std::get<7>(h);
    arma::mat A = arma::zeros(nxy, nxy); // Aerodynamic Matrix
    arma::mat L; // Lower triangular matrix
    arma::mat U; // Upper triangular matrix
    arma::mat P; // Permutation matrix
    arma::vec b = arma::zeros(nxy);
    arma::mat nC = calculateNormal(); // Unit normal vector of the wing
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
    Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi, arma::mat _x, arma::mat _y, std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec> _h)
        : chi(_chi), x(_x), y(_y), h(_h) {}
    Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi, arma::mat _x, arma::mat _y, arma::mat _z, std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec> _h)
        : chi(_chi), x(_x), y(_y), z(_z), h(_h) {}
    explicit Wing(std::array<Lagrange::CurveInterpolant*, 4> _chi) : Wing(fromTransfiniteQuadMap(_chi)) {}
    Wing(arma::mat _z, std::array<Lagrange::CurveInterpolant*, 4> _chi) : Wing(fromTransfiniteQuadMap(_z, _chi)) {}
    // Sets the dynamic pressure
    void dynamicPressure(double);
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
    void checkMesh();
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
    void output(std::string);
    void operator()(Symmetry _sym)
    {
        sym = _sym;
    }
    friend class Aerodynamics;
private:
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