#include "Airfoil.hpp"

double Airfoil::x1(double x0, double _x)
{
    return x0 - _x;
}

double Airfoil::z1(double z0, double _z)
{
    return z0 - _z;
}

double Airfoil::x2(double x0, double _x, double xi0, double _xi)
{
    return x1(x0, _x)/(xi0 - _xi);
}

double Airfoil::z2(double z0, double _z, double xi0, double _xi)
{
    return z1(z0, _z)/(xi0 - _xi);
}

double Airfoil::x3(double dx0, double x0, double _x, double xi0, double _xi)
{
    return (_x - x0 - dx0*(_xi - xi0))/pow(_xi - xi0, 2);
}

double Airfoil::z3(double dz0, double z0, double _z, double xi0, double _xi)
{
    return (_z - z0 - dz0*(_xi - xi0))/pow(_xi - xi0, 2);
}

double Airfoil::r2(double x0, double _x, double z0, double _z, double xi0, double _xi)
{
    return pow(x2(x0, _x, xi0, _xi), 2) + pow(z2(z0, _z, xi0, _xi), 2);
}

double Airfoil::k1(double dr0, double dr, double dx0, double dz0, double x0, double _x, double z0, double _z, double xi0, double _xi)
{
    return (k2(dr0, dr, dx0, dz0, x0, _x, z0, _z, xi0, _xi) + k3(dr0, dx0, dz0, x0, _x, z0, _z, xi0, _xi))/r2(x0, _x, z0, _z, xi0, _xi)/dr0;
}

double Airfoil::k2(double dr0, double dr, double dx0, double dz0, double x0, double _x, double z0, double _z, double xi0, double _xi)
{
    return (dx0*x2(x0, _x, xi0, _xi) + dz0*z2(z0, _z, xi0, _xi))*(dr - dr0)/(xi0 - _xi);
}

double Airfoil::k3(double dr0, double dx0, double dz0, double x0, double _x, double z0, double _z, double xi0, double _xi)
{
    return dr0*(x2(x0, _x, xi0, _xi)*x3(dx0, x0, _x, xi0, _xi) + z2(z0, _z, xi0, _xi)*z3(dz0, z0, _z, xi0, _xi));
}