#pragma once
#include <memory>
#include "TensorField.hpp"
#include "Metric.hpp"

class Membrane
{
    std::array<Lagrange::CurveInterpolant*, 4> chi;
    double Et = 1000; // Young's modulus multiplied by thickness
    double nu = 0.3; // Poisson's ratio
    double D = Et/(1 - nu*nu); // Material constant
    arma::mat x = arma::zeros(2, 2); // x-coordinates of nodes
    arma::mat y = arma::zeros(2, 2); // y-coordinates of nodes
    size_t nx = x.n_rows; // Number of nodes in x-direction
    size_t ny = y.n_cols; // Number of nodes in y-direction
    size_t nxy = nx*ny; // Product of nx and ny
    TensorField   z{nx, ny}; // z-coordinates of nodes
    TensorField  v1{nx, ny}; // Covariant displacement 1
    TensorField  v2{nx, ny}; // Covariant displacement 2
    TensorField n11{nx, ny}; // Contravariant stress component 11
    TensorField n12{nx, ny}; // Contravariant stress component 12
    TensorField n22{nx, ny}; // Contravariant stress component 22
    arma::mat gamma_11 = arma::zeros(nx, ny); // strain tensor (component 11)
    arma::mat gamma_12 = arma::zeros(nx, ny); // strain tensor (component 12)
    arma::mat gamma_22 = arma::zeros(nx, ny); // strain tensor (component 22)
    arma::vec x1 = Chebyshev::gaussLobatto(nx); // Cartesian x-coordinate
    arma::vec x2 = Chebyshev::gaussLobatto(ny); // Cartesian y-coordinate
    arma::mat D1  = Chebyshev::derivativeMatrix(x1, Derivative::first); // Derivative Matrix first in 1-direction
    arma::mat D2  = Chebyshev::derivativeMatrix(x2, Derivative::first); // Derivative Matrix first in 2-direction
    arma::mat D11 = Chebyshev::derivativeMatrix(x1, Derivative::second); // Derivative Matrix second in 1-direction
    arma::mat D22 = Chebyshev::derivativeMatrix(x2, Derivative::second); // Derivative Matrix second in 2-direction
    arma::field<arma::mat> J = Jacobian(x, y, D1, D2);
    arma::cube e_c = MetricCo(J); // Covariant metric tensor of the surface
    arma::cube ec  = MetricContra(e_c); // Contravariant metric tensor of the surface
    arma::mat e    = e_11()%e_22() - pow(e_12(), 2);
    arma::cube gam = Christoffel(e_c, ec, D1, D2); // Christoffel symbols
    arma::mat Psi1111 = 2*(g111()%g111() + g211()%g112()) - D1*g111();
    arma::mat Psi1112 = 2*(g111()%g211() + g211()%g212()) - D1*g211();
    arma::mat Psi2221 = 2*(g122()%g112() + g222()%g122()) - g122()*D2.t();
    arma::mat Psi2222 = 2*(g122()%g212() + g222()%g222()) - g222()*D2.t();
    arma::mat Psi1121 = 2*(g111()%g112() + g112()%g212()) - g111()*D2.t();
    arma::mat Psi1122 = 2*(g112()%g211() + g212()%g212()) - g211()*D2.t();
    arma::mat Psi2121 = 2*(g112()%g112() + g212()%g122()) - D1*g122();
    arma::mat Psi2122 = 2*(g112()%g212() + g212()%g222()) - D1*g222();
    arma::mat Psi1221 = g122()%g111() + g222()%g112() + g112()%g112() + g212()%g122() - g112()*D2.t();
    arma::mat Psi1222 = g122()%g211() + g112()%g212() + 2*g212()%g222() - g212()*D2.t();
    arma::mat Psi2111 = g111()%g212() + g211()%g222() + g112()%g211() + g212()%g212() - D1*g212();
    arma::mat Psi2112 = 2*g111()%g112() + g211()%g122() + g212()%g112() - D1*g112();
    arma::mat H1111 = e11()%e11();
    arma::mat H1112 = e11()%e12();
    arma::mat H1222 = e12()%e22();
    arma::mat H1122 = e12()%e12() + nu/e;
    arma::mat H1221 = (e11()%e22() + e12()%e12() - nu/e)/2;
    arma::mat H2222 = e22()%e22();
    arma::rowvec h_1s1_east = sqrt(e11().row(nx-1));
    arma::rowvec h_11s_east = 1/h_1s1_east;
    arma::rowvec h_22s_east = sqrt(e_22().row(nx-1));
    arma::rowvec h_2s2_east = 1/h_22s_east;
    arma::rowvec h_12s_east = e_12().row(nx-1)/sqrt(e_22().row(nx-1));
    arma::rowvec h_1s2_east = h_11s_east%e12().row(nx-1);
    arma::rowvec h_1s1_west = sqrt(e11().row(0));
    arma::rowvec h_11s_west = 1/h_1s1_west;
    arma::rowvec h_22s_west = sqrt(e_22().row(0));
    arma::rowvec h_2s2_west = 1/h_22s_west;
    arma::rowvec h_12s_west = e_12().row(0)/sqrt(e_22().row(0));
    arma::rowvec h_1s2_west = h_11s_west%e12().row(0);
    arma::vec h_11s_south = sqrt(e_11().col(0));
    arma::vec h_1s1_south = 1/h_11s_south;
    arma::vec h_2s2_south = sqrt(e22().col(0));
    arma::vec h_22s_south = 1/h_2s2_south;
    arma::vec h_21s_south = e_12().col(0)/sqrt(e_11().col(0));
    arma::vec h_2s1_south = h_22s_south%e12().col(0);
    arma::vec h_11s_north = sqrt(e_11().col(ny-1));
    arma::vec h_1s1_north = 1/h_11s_north;
    arma::vec h_2s2_north = sqrt(e22().col(ny-1));
    arma::vec h_22s_north = 1/h_2s2_north;
    arma::vec h_21s_north = e_12().col(ny-1)/sqrt(e_11().col(ny-1));
    arma::vec h_2s1_north = h_22s_north%e12().col(ny-1);
    arma::mat p = arma::zeros(nx, ny); // Pressure on the nodes
    arma::vec p1 = arma::zeros(nxy); // In-plane load in 1-direction
    arma::vec p2 = arma::zeros(nxy); // In-plane load in 2-direction
    arma::mat S = arma::zeros(nxy, nxy); // Jacobian
    arma::mat L = arma::zeros(nxy, nxy); // Lower triangular matrix
    arma::mat U = arma::zeros(nxy, nxy); // Upper triangular matrix
    arma::mat P = arma::zeros(nxy, nxy); // Permutation matrix
    arma::vec b = arma::zeros(nxy); // Right side of linear equation system
    size_t iter = 100; // Max. # of iterations for nonlinear solution (default: 1000)
    double residualTarget = 1e-10; // Target residual for nonlinear solution
    size_t substeps = 1;
    Analysis analysis = Analysis::linear;
    arma::mat DD1   = ddx();
    arma::mat DD2   = ddy();
    arma::mat D2D11 = d2dx2();
    arma::mat D2D12 = d2dxdy();
    arma::mat D2D22 = d2dy2();
    arma::mat v1__1 = arma::zeros(nx, ny);
    arma::mat v1__2 = arma::zeros(nx, ny);
    arma::mat v2__1 = arma::zeros(nx, ny);
    arma::mat v2__2 = arma::zeros(nx, ny);
    Membrane fromTransfiniteQuadMap(std::array<Lagrange::CurveInterpolant*, 4> _chi)
    {
        chi = _chi;
        auto [x, y] = Lagrange::TransfiniteQuadMap(chi);
        return {x, y};
    }
public:
    Membrane() = default;
    // Constructor to set x and y
    Membrane(arma::mat _x, arma::mat _y) : x(_x), y(_y) {}
    Membrane(std::array<Lagrange::CurveInterpolant*, 4> _chi) : Membrane(fromTransfiniteQuadMap(_chi)) {}
    // Sets the Young's modulus times thickness (Et)
    void youngsModulus(const double _Et);
    // Sets the Poisson's ratio (nu)
    void poissonsRatio(const double _nu);
    // Sets the number of iterations for nonlinear analysis (default: iter = 100)
    void iterations(size_t _iter)
    {
        iter = _iter;
    }
    // Sets the number of substeps for nonlinear analysis (default: substeps = 1)
    void substepControl(const double _substeps)
    {
        substeps = _substeps;
    }
    void checkMesh();
    // Calculates prestrain and pretension (not properly tested)
    void planeStrain();
    // Apply a constant load
    void load(const double f)
    {
        p.fill(f);
    }
    // Apply a distributed load defined by matrix
    void load(const arma::mat f)
    {
        p = f;
    }
    // Apply a distributed load defined as a function p = f(x, y)
    void load(const std::function<double(double, double)>);
    void inPlane1(const arma::mat _p1)
    {
        p1 = vectorise(_p1);
    }
    void inPlane1(const double _p1)
    {
        p1.fill(_p1);
    }
    void inPlane2(const arma::mat _p2)
    {
        p2 = vectorise(_p2);
    }
    void inPlane2(const double _p2)
    {
        p2.fill(_p2);
    }
    void inPlane11(const arma::mat _n11)
    {
        n11 = _n11;
    }
    void inPlane11(const double _n11)
    {
        n11.fill(_n11);
    }
    void inPlane12(const arma::mat _n12)
    {
        n12 = _n12;
    }
    void inPlane12(const double _n12)
    {
        n12.fill(_n12);
    }
    void inPlane22(const arma::mat _n22)
    {
        n22 = _n22;
    }
    void inPlane22(const double _n22)
    {
        n22.fill(_n22);
    }
    void inPlaneTensorial(const arma::vec, const arma::vec);
    void inPlaneKartesian(const arma::vec, const arma::vec);
    void inPlaneKartesian(const arma::mat, const arma::mat, const arma::mat);
    arma::mat X() const
    {
        return x;
    }
    double X(const size_t i, const size_t j) const
    {
        return x(i, j);
    }
    arma::mat Y() const
    {
        return y;
    }
    double Y(const size_t i, const size_t j) const
    {
        return y(i, j);
    }
    arma::mat Z() const
    {
        return z;
    }
    double Z(const size_t i, const size_t j) const
    {
        return z(i, j);
    }
    arma::mat V1() const
    {
        return v1;
    }
    double V1(const size_t i, const size_t j) const
    {
        return v1(i, j);
    }
    arma::mat V2() const
    {
        return v2;
    }
    double V2(const size_t i, const size_t j) const
    {
        return v2(i, j);
    }
    arma::mat N11() const
    {
        return n11;
    }
    arma::mat N12() const
    {
        return n12;
    }
    arma::mat N22() const
    {
        return n22;
    }
    arma::mat J11() const
    {
        return J(0, 0);
    }
    double J11(const size_t i, const size_t j) const
    {
        return J(0, 0)(i, j);
    }
    arma::mat J12() const
    {
        return J(0, 1);
    }
    double J12(const size_t i, const size_t j) const
    {
        return J(0, 1)(i, j);
    }
    arma::mat J21() const
    {
        return J(1, 0);
    }
    double J21(const size_t i, const size_t j) const
    {
        return J(1, 0)(i, j);
    }
    arma::mat J22() const
    {
        return J(1, 1);
    }
    double J22(const size_t i, const size_t j) const
    {
        return J(1, 1)(i, j);
    }
    arma::mat e_11() const
    {
        return e_c.slice(0);
    }
    double e_11(const size_t i, const size_t j) const
    {
        return e_c(i, j, 0);
    }
    arma::mat e_12() const
    {
        return e_c.slice(1);
    }
    double e_12(const size_t i, const size_t j) const
    {
        return e_c(i, j, 1);
    }
    arma::mat e_22() const
    {
        return e_c.slice(2);
    }
    double e_22(const size_t i, const size_t j) const
    {
        return e_c(i, j, 2);
    }
    arma::mat e11() const
    {
        return ec.slice(0);
    }
    double e11(const size_t i, const size_t j) const
    {
        return ec(i, j, 0);
    }
    arma::mat e12() const
    {
        return ec.slice(1);
    }
    double e12(const size_t i, const size_t j) const
    {
        return ec(i, j, 1);
    }
    arma::mat e22() const
    {
        return ec.slice(2);
    }
    double e22(const size_t i, const size_t j) const
    {
        return ec(i, j, 2);
    }
    arma::mat g111() const
    {
        return gam.slice(0);
    }
    double g111(const size_t i, const size_t j) const
    {
        return gam(i, j, 0);
    }
    arma::mat g112() const
    {
        return gam.slice(1);
    }
    double g112(const size_t i, const size_t j) const
    {
        return gam(i, j, 1);
    }
    arma::mat g122() const
    {
        return gam.slice(2);
    }
    double g122(const size_t i, const size_t j) const
    {
        return gam(i, j, 2);
    }
    arma::mat g211() const
    {
        return gam.slice(3);
    }
    double g211(const size_t i, const size_t j) const
    {
        return gam(i, j, 3);
    }
    arma::mat g212() const
    {
        return gam.slice(4);
    }
    double g212(const size_t i, const size_t j) const
    {
        return gam(i, j, 4);
    }
    arma::mat g222() const
    {
        return gam.slice(5);
    }
    double g222(const size_t i, const size_t j) const
    {
        return gam(i, j, 5);
    }
    // Set homogeneous boundary condition at specified boundary with specified type
    void boundary(const Field field, const Direction dir, const BC bc)
    {
        boundary(field, dir, bc, 0);
    }
    void boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc)
    {
        boundary(field, dir, bc, 0);
    }
    // Set arbitrary boundary condition at specified boundary with specified type
    template <class C> void boundary(const Field, const Direction, const BC, const C);
    template <class C> void boundary(const Field, const Lagrange::CurveInterpolant*, const BC, const C);
    // Set boundary condition with constants for homogeneous boundary condition
    void boundary(const Field field, const Direction dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    void boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2)
    {
        boundary(field, dir, bc, _r1, _r2, 0);
    }
    // Set boundary condition with constants for Robin boundary condition
    template <class C> void boundary(const Field, const Direction, const BC, const double, const double, const C);
    template <class C> void boundary(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const C);
    // Linear analysis
    void linear();
    // Semilinear analysis
    void semilinear();
    // Nonlinear analysis
    void nonlinear();
    // Calculates the integral of a given matrix over the surface of the membrane
    double integrate(arma::mat);
    // Calculates the integral of a given field over the surface of the membrane
    double integrate(const Field);
    // Calculates the elastic potential or strain energy per unit area of the middle surface
    double elasticPotential();
    std::pair<arma::mat, arma::mat> kartesianDisplacements();
    void principalStresses(std::string, std::string);
    void principalStrains(std::string, std::string, std::string);
    // Output y, x and a chosen field to chosen file
    void output(const std::string, const Field);
    double operator()(const size_t, const size_t, const Field);
    double operator()(const Field field, const size_t i, const size_t j)
    {
        return operator()(i, j, field);
    }
    void operator()(const arma::mat, const arma::mat);
    void operator()(const std::array<Lagrange::CurveInterpolant*, 4>);
    TensorField operator()(const Field);
    friend class Structure;
private:
    arma::mat ddx()
    {
        return arma::kron(arma::eye(ny, ny), D1);
    }
    arma::mat ddx(const arma::mat);
    arma::mat ddy()
    {
        return arma::kron(D2, arma::eye(nx, nx));
    }
    arma::mat ddy(const arma::mat);
    arma::mat d2dx2()
    {
        return arma::kron(arma::eye(ny, ny), D11);
    }
    arma::mat d2dx2(const arma::mat);
    arma::mat d2dy2()
    {
        return arma::kron(D22, arma::eye(nx, nx));
    }
    arma::mat d2dy2(const arma::mat);
    arma::mat d2dxdy()
    {
        return arma::kron(D2, D1);
    }
    arma::mat d2dxdy(const arma::mat);
    arma::mat constant(const arma::mat H)
    {
        return arma::diagmat(arma::vectorise(H));
    }
    std::unique_ptr<TensorField> setField(const Field);
    void structuralMatrix();
    void solve_S();
    void solve_b();
    double armijoSemilinear(arma::vec, arma::vec&);
    double armijoNonlinear(arma::vec, arma::mat&, arma::mat&, arma::mat&);
    double residualLevelFunctionSemilinear(arma::vec, arma::vec&);
    double residualLevelFunctionNonlinear(arma::vec, arma::mat&, arma::mat&, arma::mat&);
    void zBoundary(const BC, const double, const size_t, const size_t,
                   const double, const double, const double, const double);
    void v1BoundaryWestEast(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                            const double, const double, const double, const double);
    void v2BoundaryWestEast(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                            const double, const double, const double, const double, const double);
    void v1BoundarySouthNorth(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                              const double, const double, const double, const double, const double);
    void v2BoundarySouthNorth(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                              const double, const double, const double, const double);
    void n11BoundaryLinear(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&, const double);
    void n12BoundaryLinear(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                           const double, const double, const double);
    void n21BoundaryLinear(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                           const double, const double, const double);
    void n22BoundaryLinear(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&, const double);
    void zBoundarySemilinear(const BC, const double, const size_t, const size_t, arma::mat&, arma::vec&,
                             const double, const double, const double, const double, const double, const double);
    void zBoundaryNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                             const double, const double, const double, const double, const double, const double);
    void v1BoundaryWestEastNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                                     const double, const double, const double, const double);
    void v2BoundaryWestEastNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                                     const double, const double, const double, const double, const double);
    void v1BoundarySouthNorthNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                                       const double, const double, const double, const double, const double);
    void v2BoundarySouthNorthNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                                       const double, const double, const double, const double);
    void n11BoundaryNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                              const double, const double, const double);
    void n12BoundaryNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                              const double, const double, const double, const double, const double);
    void n21BoundaryNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                              const double, const double, const double, const double, const double);
    void n22BoundaryNonlinear(const BC, const double, const size_t, const size_t, arma::mat&,
                              const double, const double, const double);
};