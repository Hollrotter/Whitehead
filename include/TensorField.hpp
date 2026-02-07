#pragma once
#include "enums.hpp"
#include "misc.hpp"

class TensorField : public arma::mat
{
    BC  westBC = BC::None; // Boundary type of wester boundary
    BC  eastBC = BC::None; // Boundary type of eastern boundary
    BC northBC = BC::None; // Boundary type of northern boundary
    BC southBC = BC::None; // Boundary type of southern boundary
    arma::vec  west; // Boundary value for western boundary
    arma::vec  east; // Boundary value for eastern boundary
    arma::vec north; // Boundary value for northern boundary
    arma::vec south; // Boundary value for southern boundary
    double r1West; // First constant for Robin boundary condition
    double r2West; // Second constant for Robin boundary condition
    double r1East; // First constant for Robin boundary condition
    double r2East; // Second constant for Robin boundary condition
    double r1North; // First constant for Robin boundary condition
    double r2North; // Second constant for Robin boundary condition
    double r1South; // First constant for Robin boundary condition
    double r2South; // Second constant for Robin boundary condition
public:
    TensorField(size_t n, size_t m) : west(arma::ones(m)), east(arma::ones(m)), north(arma::ones(n)), south(arma::ones(n))
    {
        zeros(n, m);
    }
    void initBC(size_t n, size_t m)
    {
        west  = arma::ones(m);
        east  = arma::ones(m);
        north = arma::ones(n);
        south = arma::ones(n);
    }
    TensorField& operator=(arma::mat A)
    {
        arma::mat::operator=(A);
        return *this;
    }
    inline arma::mat &operator+=(const arma::mat &A)
    {
        return arma::mat::operator+=(A);
    }
    inline arma::mat &operator-=(const arma::mat &A)
    {
        return arma::mat::operator-=(A);
    }
    template <Number C> arma::mat &operator*=(const C &c)
    {
        return arma::mat::operator*=(c);
    }
    inline arma::mat &operator*=(const arma::mat &A)
    {
        return arma::mat::operator*=(A);
    }
    inline arma::mat &operator/=(const arma::mat &A)
    {
        return arma::mat::operator/=(A);
    }
    inline arma::mat &operator/=(const double &d)
    {
        return arma::mat::operator/=(d);
    }
    friend class Membrane;
};

inline arma::mat  operator+(TensorField A, const TensorField &B)
{
    return A += B;
}
inline arma::mat  operator-(TensorField A, const TensorField &B)
{
    return A -= B;
}
template <Number C> arma::mat operator*(C c, const TensorField &T)
{
    return static_cast<arma::mat>(T) *= c;
}
template <Number C> arma::mat operator*(const TensorField &T, C c)
{
    return static_cast<arma::mat>(T) *= c;
}
inline arma::mat operator*(arma::mat M, TensorField T)
{
    return M *= static_cast<arma::mat>(T);
}
inline arma::mat operator*(TensorField T, arma::mat M)
{
    return static_cast<arma::mat>(T) *= M;
}
inline arma::mat operator%(arma::mat A, TensorField T)
{
    return A % static_cast<arma::mat>(T);
}
inline arma::mat operator%(TensorField T, arma::mat A)
{
    return static_cast<arma::mat>(T) % A;
}
inline arma::mat operator%(TensorField T1, TensorField T2)
{
    return static_cast<arma::mat>(T1) % static_cast<arma::mat>(T2);
}
inline arma::mat operator/(TensorField T, arma::mat M)
{
    return static_cast<arma::mat>(T) /= M;
}
inline arma::mat operator/(TensorField T, double d)
{
    return static_cast<arma::mat>(T) /= d;
}
inline arma::mat operator/(arma::mat A, TensorField T)
{
    return A /= static_cast<arma::mat>(T);
}