#include "Membrane.hpp"

arma::mat Membrane::ddx(arma::mat H)
{
    arma::mat DDX = arma::kron(arma::eye(ny, ny), D1);
    for (size_t j = 0, k = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++, k++)
            DDX.row(k) *= H(i, j);
    return DDX;
}

arma::mat Membrane::ddy(arma::mat H)
{
    arma::mat DDY = arma::kron(D2, arma::eye(nx, nx));
    for (size_t j = 0, k = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++, k++)
            DDY.row(k) *= H(i, j);
    return DDY;
}

arma::mat Membrane::d2dx2(arma::mat H)
{
    arma::mat D2DX2 = arma::kron(arma::eye(ny, ny), D11);
    for (size_t j = 0, k = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++, k++)
            D2DX2.row(k) *= H(i, j);
    return D2DX2;
}

arma::mat Membrane::d2dy2(arma::mat H)
{
    arma::mat D2DY2 = arma::kron(D22, arma::eye(nx, nx));
    for (size_t j = 0, k = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++, k++)
            D2DY2.row(k) *= H(i, j);
    return D2DY2;
}

arma::mat Membrane::d2dxdy(arma::mat H)
{
    arma::mat D2DXDY = arma::kron(D2, D1);
    for (size_t j = 0, k = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++, k++)
            D2DXDY.row(k) *= H(i, j);
    return D2DXDY;
}