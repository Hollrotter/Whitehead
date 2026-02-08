#pragma once
#include "Chebyshev.hpp"
#include "Point.hpp"

class Membrane;
class Structure;

namespace Lagrange
{
    // Calculates the barycentric weights needed for the Lagrange-Interpolation
    arma::vec barycentricWeights(const arma::vec x);
    arma::vec barycentricWeights(const size_t n);
    double interpolation(const double x, const arma::vec X, const arma::vec f, const arma::vec w);
    arma::vec interpolation(const arma::vec x, const arma::vec X, const arma::vec f, const arma::vec w);
    double interpolantDerivative(const double x, const arma::vec X, const arma::vec f, const arma::vec w);
    arma::vec interpolantDerivative(const arma::vec x, const arma::vec X, const arma::vec f, const arma::vec w);
    // Calculates the Matrix for the Lagrange-Interpolation
    arma::mat interpolationMatrix(const arma::vec x, const arma::vec xi);
    // Calculates the Matrix for the Lagrange-Interpolation for x = f(x1, x2) and y = f(x1, x2)
    arma::mat interpolationMatrix(const arma::vec x1, const arma::vec x2, const arma::mat x, const arma::mat y);
    // Performs the Lagrange-Interpolation in 2D given the Matrices for x and y directions.
    arma::mat interpolation2D(const arma::mat Tx, const arma::mat Ty, const arma::mat z, const arma::vec xi, const arma::vec eta);
    class CurveInterpolant
    {
        arma::vec x = arma::zeros(1);
        arma::vec y = arma::zeros(1);
        arma::vec nodes;
        arma::vec w = barycentricWeights(nodes);
        CurveType curveType = CurveType::Boundary;
        CurveInterpolant(arma::vec _x, arma::vec _y, double r) : x(_x), y(_y), nodes(Chebyshev::gaussLobatto(x.size())) {}
        CurveInterpolant arc(Point p1, Point p2, Point pm, double r, size_t n)
        {
            double phi1 = atan(p1.Y()/p1.X());
            double phi2 = atan(p2.Y()/p2.X());
            arma::vec phi = phi1 + (phi2-phi1)*(1+Chebyshev::gaussLobatto(n))/2;
            arma::vec x = r*cos(phi) + pm.X();
            arma::vec y = r*sin(phi) + pm.Y();
            return {x, y, r};
        }
    public:
        CurveInterpolant() = default;
        CurveInterpolant(arma::vec _x, arma::vec _y) : x(_x), y(_y), nodes(parametrize()) {}
        CurveInterpolant(Point p1, Point p2, size_t n) : x((p2.X()-p1.X())*(Chebyshev::gaussLobatto(n)+1)/2 + p1.X()),
                                                         y((p2.Y()-p1.Y())*(Chebyshev::gaussLobatto(n)+1)/2 + p1.Y()), nodes(Chebyshev::gaussLobatto(n)) {};
        CurveInterpolant(Point p1, Point p2, Point pm, double r, size_t n) : CurveInterpolant(arc(p1, p2, pm, r, n)) {};
        template <class S> std::pair<S, S> evaluate(const S s)
        {
            return {interpolation(s, nodes, x, w), interpolation(s, nodes, y, w)};
        }
        template <class S> std::pair<S, S> derivative(const S s)
        {
            return {interpolantDerivative(s, nodes, x, w), interpolantDerivative(s, nodes, y, w)};
        }
        double arclengthAt(const double s, const arma::vec X)
        {
            return interpolation(s, nodes, X, w);
        }
        double operator()(const size_t i)
        {
            return nodes(i);
        }
        arma::vec getNodes()
        {
            return nodes;
        }
        friend class ::Membrane;
        friend class ::Structure;
    private:
        arma::vec parametrize();
    };
    std::pair<arma::mat, arma::mat> TransfiniteQuadMap(const std::array<CurveInterpolant*, 4> chi);
    std::tuple<arma::mat, arma::mat, arma::mat, arma::mat> TransfiniteQuadMetrics(const std::array<CurveInterpolant*, 4> chi);
}