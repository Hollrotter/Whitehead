#include "Membrane.hpp"

/**
 * @brief 
 * 
 * @tparam C 
 * @param field 
 * @param dir 
 * @param bc 
 * @param val 
 */
template <class C> void Membrane::boundary(const Field field, const Direction dir, const BC bc, const C val)
{
    std::unique_ptr<TensorField> f = setField(field);
    switch(dir)
    {
        case Direction::N:
            f->northBC = bc;
            f->north   = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->north*val);
            break;
        case Direction::S:
            f->southBC = bc;
            f->south   = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->south*val);
            break;
        case Direction::W:
            f->westBC  = bc;
            f->west    = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->west*val);
            break;
        case Direction::E:
            f->eastBC  = bc;
            f->east    = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->east*val);
            break;
    }
    f.release();
}

template <class C> void Membrane::boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const C val)
{
    std::unique_ptr<TensorField> f = setField(field);
    if(dir == chi[0])
    {
        f->southBC = bc;
        f->south   = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->south*val);
    }
    else if(dir == chi[1])
    {
        f->eastBC  = bc;
        f->east    = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->east*val);
    }
    else if(dir == chi[2])
    {
        f->northBC = bc;
        f->north   = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->north*val);
    }
    else if(dir == chi[3])
    {
        f->westBC  = bc;
        f->west    = (typeid(val) == typeid(arma::vec) ? arma::vec(val) : f->west*val);
    }
    else
    {
        std::println("Not a curve of the membrane!");
        f.release();
        exit(EXIT_FAILURE);
    }
    f.release();
}

/**
 * @brief 
 * 
 * @tparam C 
 * @param field 
 * @param dir 
 * @param bc 
 * @param _r1 
 * @param _r2 
 * @param val 
 */
template <class C> void Membrane::boundary(const Field field, const Direction dir, const BC bc, const double _r1, const double _r2, const C val)
{
    std::unique_ptr<TensorField> f = setField(field);
    switch(dir)
    {
        case Direction::N:
            f->r1North = _r1;
            f->r2North = _r2;
            break;
        case Direction::S:
            f->r1South = _r1;
            f->r2South = _r2;
            break;
        case Direction::W:
            f->r1West = _r1;
            f->r2West = _r2;
            break;
        case Direction::E:
            f->r1East = _r1;
            f->r2East = _r2;
            break;
    }
    f.release();
    boundary(field, dir, bc, val);
}

template <class C> void Membrane::boundary(const Field field, const Lagrange::CurveInterpolant* dir, const BC bc, const double _r1, const double _r2, const C val)
{
    std::unique_ptr<TensorField> f = setField(field);
    if(dir == chi[0])
    {
        f->r1South = _r1;
        f->r2South = _r2;
    }
    else if(dir == chi[1])
    {
        f->r1East = _r1;
        f->r2East = _r2;
    }
    else if(dir == chi[2])
    {
        f->r1North = _r1;
        f->r2North = _r2;
    }
    else if(dir == chi[3])
    {
        f->r1West = _r1;
        f->r2West = _r2;
    }
    else
    {
        std::println("Not a curve of the membrane!");
        f.release();
        exit(EXIT_FAILURE);
    }
    f.release();
    boundary(field, dir, bc, val);
}

void Membrane::zBoundary(const BC bc, const double val, const size_t i, const size_t j,
                         const double r1, const double r2, const double h1, const double h2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            S.row(k).zeros();
            S(k, k) = 1;
            b(k) = val;
            break;
        case BC::Neumann:
            S.row(k) = DD1.row(k)*h1 + DD2.row(k)*h2;
            b(k) = val;
            break;
        case BC::Robin:
            S.row(k) = r2*(DD1.row(k)*h1 + DD2.row(k)*h2);
            S(k, k) += r1;
            b(k) = val;
            break;
        case BC::None:
            break;
    }
}

void Membrane::v1BoundaryWestEast(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                  const double r1, const double r2, const double h_1s1, const double h_1s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            bv(k) = val;
            A.row(k).zeros();
            A(k, k)     = h_1s1;
            A(k, k+nxy) = h_1s2;
            break;
        case BC::Neumann:
        {
            bv(k) = val;
            double c11 = pow(h_1s1, 2);
            double c12 = h_1s1*h_1s2;
            double c22 = pow(h_1s2, 2);
            A.row(k).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k, k)     -= c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2);
            A(k, k+nxy) -= c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5);
            break;
        }
        case BC::Robin:
        {
            bv(k) = val;
            double c11 = pow(h_1s1, 2);
            double c12 = h_1s1*h_1s2;
            double c22 = pow(h_1s2, 2);
            A.row(k).cols(  0,   nxy-1) = r2*(c11*DD1.row(k) + c12*DD2.row(k));
            A.row(k).cols(nxy, 2*nxy-1) = r2*(c12*DD1.row(k) + c22*DD2.row(k));
            A(k, k)     += r1*h_1s1 - r2*(c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2));
            A(k, k+nxy) += r1*h_1s2 - r2*(c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::v2BoundaryWestEast(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                    const double r1, const double r2, const double h_1s1, const double h_1s2, const double h_2s2)
{
    size_t k = i + j*nx;
        switch (bc)
        {
            case BC::Dirichlet:
                bv(k+nxy) = val;
                A.row(k+nxy).zeros();
                A(k+nxy, k+nxy) = h_2s2;
                break;
            case BC::Neumann:
            {
                bv(k+nxy) = val;
                double c1 = h_2s2*h_1s1;
                double c2 = h_2s2*h_1s2;
                A.row(k+nxy).zeros();
                A.row(k+nxy).cols(nxy, 2*nxy-1) = c1*DD1.row(k) + c2*DD2.row(k);
                A(k+nxy, k)     -= c1*gam(i, j, 1) + c2*gam(i, j, 2);
                A(k+nxy, k+nxy) -= c1*gam(i, j, 4) + c2*gam(i, j, 5);
                break;
            }
            case BC::Robin:
            {
                bv(k+nxy) = val;
                double c1 = h_2s2*h_1s1;
                double c2 = h_2s2*h_1s2;
                A.row(k+nxy).zeros();
                A.row(k+nxy).cols(nxy, 2*nxy-1) = r2*(c1*DD1.row(k) + c2*DD2.row(k));
                A(k+nxy, k)     -=            r2*(c1*gam(i, j, 1) + c2*gam(i, j, 2));
                A(k+nxy, k+nxy) += r1*h_2s2 - r2*(c1*gam(i, j, 4) + c2*gam(i, j, 5));
                break;
            }
            case BC::None:
                break;
        }
}

void Membrane::v1BoundarySouthNorth(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                    const double r1, const double r2, const double h_1s1, const double h_2s1, const double h_2s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            bv(k) = val;
            A.row(k).zeros();
            A(k, k) = h_1s1;
            break;
        case BC::Neumann:
        {
            bv(k) = val;
            double c1 = h_1s1*h_2s1;
            double c2 = h_1s1*h_2s2;
            A.row(k).zeros();
            A.row(k).cols(0, nxy-1) = c1*DD1.row(k) + c2*DD2.row(k);
            A(k, k)     -= c1*gam(i, j, 0) + c2*gam(i, j, 1);
            A(k, k+nxy) -= c1*gam(i, j, 3) + c2*gam(i, j, 4);
            break;
        }
        case BC::Robin:
        {
            bv(k) = val;
            double c1 = h_1s1*h_2s1;
            double c2 = h_1s1*h_2s2;
            A.row(k).zeros();
            A.row(k).cols(0, nxy-1) = r2*(c1*DD1.row(k) + c2*DD2.row(k));
            A(k, k)     += r1*h_1s1 - r2*(c1*gam(i, j, 0) + c2*gam(i, j, 1));
            A(k, k+nxy) -=            r2*(c1*gam(i, j, 3) + c2*gam(i, j, 4));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::v2BoundarySouthNorth(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                    const double r1, const double r2, const double h_2s1, const double h_2s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            bv(k+nxy) = val;
            A.row(k+nxy).zeros();
            A(k+nxy, k)     = h_2s1;
            A(k+nxy, k+nxy) = h_2s2;
            break;
        case BC::Neumann:
        {
            bv(k+nxy) = val;
            double c11 = pow(h_2s1, 2);
            double c12 = h_2s1*h_2s2;
            double c22 = pow(h_2s2, 2);
            A.row(k+nxy).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k+nxy).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k+nxy, k)     -= c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2);
            A(k+nxy, k+nxy) -= c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5);
            break;
        }
        case BC::Robin:
        {
            bv(k+nxy) = val;
            double c11 = pow(h_2s1, 2);
            double c12 = h_2s1*h_2s2;
            double c22 = pow(h_2s2, 2);
            A.row(k+nxy).cols(  0,   nxy-1) = r2*(c11*DD1.row(k) + c12*DD2.row(k));
            A.row(k+nxy).cols(nxy, 2*nxy-1) = r2*(c12*DD1.row(k) + c22*DD2.row(k));
            A(k+nxy, k)     += r1*h_2s1 - r2*(c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2));
            A(k+nxy, k+nxy) += r1*h_2s2 - r2*(c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::n11BoundaryLinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv, const double h_11s)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            bv(k) = val/D/pow(h_11s, 2);
            A.row(k).cols(  0,   nxy-1) = H1111(i, j)*DD1.row(k) + H1112(i, j)*DD2.row(k);
            A.row(k).cols(nxy, 2*nxy-1) = H1112(i, j)*DD1.row(k) + H1122(i, j)*DD2.row(k);
            A(k, k)     -= gam(i, j, 0)*H1111(i, j) + 2*gam(i, j, 1)*H1112(i, j) + gam(i, j, 2)*H1122(i, j);
            A(k, k+nxy) -= gam(i, j, 3)*H1111(i, j) + 2*gam(i, j, 4)*H1112(i, j) + gam(i, j, 5)*H1122(i, j);
            break;
        }
        case BC::Neumann:
            std::println("n11 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::Robin:
            std::println("n11 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n12BoundaryLinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                 const double h_11s, const double h_12s, const double h_22s)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            bv(k+nxy) = val/D;
            double c11 = h_11s*h_12s*H1111(i, j) + h_11s*h_22s*H1112(i, j);
            double c12 = h_11s*h_12s*H1112(i, j) + h_11s*h_22s*H1221(i, j);
            double c22 = h_11s*h_12s*H1122(i, j) + h_11s*h_22s*H1222(i, j);
            A.row(k+nxy).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k+nxy).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k+nxy, k)     -= gam(i, j, 0)*c11 + 2*gam(i, j, 1)*c12 + gam(i, j, 2)*c22;
            A(k+nxy, k+nxy) -= gam(i, j, 3)*c11 + 2*gam(i, j, 4)*c12 + gam(i, j, 5)*c22;
            break;
        }
        case BC::Neumann:
            std::println("n12 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::Robin:
            std::println("n12 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n21BoundaryLinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv,
                                 const double h_11s, const double h_21s, const double h_22s)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            bv(k) = val/D;
            double c11 = h_11s*h_22s*H1112(i, j) + h_21s*h_22s*H1122(i, j);
            double c12 = h_11s*h_22s*H1221(i, j) + h_21s*h_22s*H1222(i, j);
            double c22 = h_11s*h_22s*H1222(i, j) + h_21s*h_22s*H2222(i, j);
            A.row(k).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k, k)     -= gam(i, j, 0)*c11 + 2*gam(i, j, 1)*c12 + gam(i, j, 2)*c22;
            A(k, k+nxy) -= gam(i, j, 3)*c11 + 2*gam(i, j, 4)*c12 + gam(i, j, 5)*c22;
            break;
        }
        case BC::Neumann:
            std::println("n21 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
            break;
        case BC::Robin:
            std::printf("n21 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n22BoundaryLinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &bv, const double h_22s)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            bv(k+nxy) = val/D/pow(h_22s, 2);
            A.row(k+nxy).cols(  0,   nxy-1) = H1122(i, j)*DD1.row(k) + H1222(i, j)*DD2.row(k);
            A.row(k+nxy).cols(nxy, 2*nxy-1) = H1222(i, j)*DD1.row(k) + H2222(i, j)*DD2.row(k);
            A(k+nxy, k)     -= gam(i, j, 0)*H1122(i, j) + 2*gam(i, j, 1)*H1222(i, j) + gam(i, j, 2)*H2222(i, j);
            A(k+nxy, k+nxy) -= gam(i, j, 3)*H1122(i, j) + 2*gam(i, j, 4)*H1222(i, j) + gam(i, j, 5)*H2222(i, j);
            break;
        }
        case BC::Neumann:
            std::println("n22 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
            break;
        case BC::Robin:
            std::println("n22 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::zBoundarySemilinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A, arma::vec &b,
                                     const double r1, const double r2, const double z_1, const double z_2, const double h1, const double h2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k) = z(i, j) - val;
            A.row(k).zeros();
            A(k, k) = 1;
            break;
        case BC::Neumann:
            b(k) = z_1*h1 + z_2*h2 - val;
            A.row(k) = DD1.row(k)*h1 + DD2.row(k)*h2;
            break;
        case BC::Robin:
            b(k) = r1*z(i, j) + r2*(z_1*h1 + z_2*h2) - val;
            A.row(k) = r2*(DD1.row(k)*h1 + DD2.row(k)*h2);
            A(k, k) += r1;
            break;
        case BC::None:
            break;
    }
}

void Membrane::zBoundaryNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                  const double r1, const double r2, const double z_1, const double z_2, const double h1, const double h2)
{
    size_t k = i + j*nx;
    size_t nxy2 = 2*nxy;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k+nxy2) = z(i, j) - val;
            A.row(k+nxy2).zeros();
            A(k+nxy2, k+nxy2) = 1;
            break;
        case BC::Neumann:
            b(k+nxy2) = z_1*h1 + z_2*h2 - val;
            A.row(k+nxy2).zeros();
            A.row(k+nxy2).cols(nxy2, 3*nxy-1) = DD1.row(k)*h1 + DD2.row(k)*h2;
            break;
        case BC::Robin:
            b(k+nxy2) = r1*z(i, j) + r2*(z_1*h1 + z_2*h2) - val;
            A.row(k+nxy2).zeros();
            A.row(k+nxy2).cols(nxy2, 3*nxy-1) = r2*(DD1.row(k)*h1 + DD2.row(k)*h2);
            A(k+nxy2, k+nxy2) += r1;
            break;
        case BC::None:
            break;
    }
}

void Membrane::v1BoundaryWestEastNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                           const double r1, const double r2, const double h_1s1, const double h_1s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k) = h_1s1*v1(i, j) + h_1s2*v2(i, j) - val;
            A.row(k).zeros();
            A(k, k)     = h_1s1;
            A(k, k+nxy) = h_1s2;
            break;
        case BC::Neumann:
        {
            double c11 = pow(h_1s1, 2);
            double c12 = h_1s1*h_1s2;
            double c22 = pow(h_1s2, 2);
            b(k) = c11*v1__1(i, j) + c12*(v1__2(i, j) + v2__1(i, j)) + c22*v2__2(i, j) - val;
            A.row(k).zeros();
            A.row(k).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k, k)     -= c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2);
            A(k, k+nxy) -= c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5);
            break;
        }
        case BC::Robin:
        {
            double c11 = pow(h_1s1, 2);
            double c12 = h_1s1*h_1s2;
            double c22 = pow(h_1s2, 2);
            b(k) = r1*(h_1s1*v1(i, j) + h_1s2*v2(i, j)) + r2*(c11*v1__1(i, j) + c12*(v1__2(i, j) + v2__1(i, j)) + c22*v2__2(i, j)) - val;
            A.row(k).zeros();
            A.row(k).cols(  0,   nxy-1) = r2*(c11*DD1.row(k) + c12*DD2.row(k));
            A.row(k).cols(nxy, 2*nxy-1) = r2*(c12*DD1.row(k) + c22*DD2.row(k));
            A(k, k)     += r1*h_1s1 - r2*(c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2));
            A(k, k+nxy) += r1*h_1s2 - r2*(c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::v2BoundaryWestEastNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                             const double r1, const double r2, const double h_1s1, const double h_1s2, const double h_2s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k+nxy) = h_2s2*v2(i, j) - val;
            A.row(k+nxy).zeros();
            A(k+nxy, k+nxy) = h_2s2;
            break;
        case BC::Neumann:
        {
            double c1 = h_1s1*h_2s2;
            double c2 = h_1s2*h_2s2;
            b(k+nxy) = c1*v2__1(i, j) + c2*v2__2(i, j) - val;
            A.row(k+nxy).zeros();
            A.row(k+nxy).cols(nxy, 2*nxy-1) = c1*DD1.row(k) + c2*DD2.row(k);
            A(k+nxy, k)     -= c1*gam(i, j, 1) + c2*gam(i, j, 2);
            A(k+nxy, k+nxy) -= c1*gam(i, j, 4) + c2*gam(i, j, 5);
            break;
        }
        case BC::Robin:
        {
            double c1 = h_1s1*h_2s2;
            double c2 = h_1s2*h_2s2;
            b(k+nxy) = r1*h_2s2*v2(i, j) + r2*(c1*v2__1(i, j) + c2*v2__2(i, j)) - val;
            A.row(k+nxy).zeros();
            A.row(k+nxy).cols(nxy, 2*nxy-1) = r2*(c1*DD1.row(k) + c2*DD2.row(k));
            A(k+nxy, k)     -=            r2*(c1*gam(i, j, 1) + c2*gam(i, j, 2));
            A(k+nxy, k+nxy) += r1*h_2s2 - r2*(c1*gam(i, j, 4) + c2*gam(i, j, 5));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::v1BoundarySouthNorthNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                           const double r1, const double r2, const double h_1s1, const double h_2s1, const double h_2s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k) = h_1s1*v1(i, j) - val;
            A.row(k).zeros();
            A(k, k) = h_1s1;
            break;
        case BC::Neumann:
        {
            double c1 = h_1s1*h_2s1;
            double c2 = h_1s1*h_2s2;
            b(k) = c1*v1__1(i, j) + c2*v1__2(i, j) - val;
            A.row(k).zeros();
            A.row(k).cols(0, nxy-1) = c1*DD1.row(k) + c2*DD2.row(k);
            A(k, k)     -= c1*gam(i, j, 0) + c2*gam(i, j, 1);
            A(k, k+nxy) -= c1*gam(i, j, 3) + c2*gam(i, j, 4);
            break;
        }
        case BC::Robin:
        {
            double c1 = h_1s1*h_2s1;
            double c2 = h_1s1*h_2s2;
            b(k) = r1*h_1s1*v1(i, j)
                 + r2*(c1*v1__1(i, j) + c2*v1__2(i, j)) - val;
            A.row(k).zeros();
            A.row(k).cols(0, nxy-1) = r2*(c1*DD1.row(k) + c2*DD2.row(k));
            A(k, k)     += r1*h_1s1 - r2*(c1*gam(i, j, 0) + c2*gam(i, j, 1));
            A(k, k+nxy) -=            r2*(c1*gam(i, j, 3) + c2*gam(i, j, 4));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::v2BoundarySouthNorthNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                             const double r1, const double r2, const double h_2s1, const double h_2s2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
            b(k+nxy) = h_2s1*v1(i, j) + h_2s2*v2(i, j) - val;
            A.row(k+nxy).zeros();
            A(k+nxy, k)     = h_2s1;
            A(k+nxy, k+nxy) = h_2s2;
            break;
        case BC::Neumann:
        {
            double c11 = pow(h_2s1, 2);
            double c12 = h_2s1*h_2s2;
            double c22 = pow(h_2s2, 2);
            b(k+nxy) = c11*v1__1(i, j) + c12*(v1__2(i, j) + v2__1(i, j)) + c22*v2__2(i, j) - val;
            A.row(k+nxy).zeros();
            A.row(k+nxy).cols(  0,   nxy-1) = c11*DD1.row(k) + c12*DD2.row(k);
            A.row(k+nxy).cols(nxy, 2*nxy-1) = c12*DD1.row(k) + c22*DD2.row(k);
            A(k+nxy, k)     -= c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2);
            A(k+nxy, k+nxy) -= c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5);
            break;
        }
        case BC::Robin:
        {
            double c11 = pow(h_2s1, 2);
            double c12 = h_2s1*h_2s2;
            double c22 = pow(h_2s2, 2);
            b(k+nxy) = r1*(h_2s1*v1(i, j) + h_2s2*v2(i, j)) + r2*(c11*v1__1(i, j) + c12*(v1__2(i, j) + v2__1(i, j)) + c22*v2__2(i, j)) - val;
            A.row(k+nxy).zeros();
            A.row(k+nxy).cols(  0,   nxy-1) = r2*(c11*DD1.row(k) + c12*DD2.row(k));
            A.row(k+nxy).cols(nxy, 2*nxy-1) = r2*(c12*DD1.row(k) + c22*DD2.row(k));
            A(k+nxy, k)     += r1*h_2s1 - r2*(c11*gam(i, j, 0) + 2*c12*gam(i, j, 1) + c22*gam(i, j, 2));
            A(k+nxy, k+nxy) += r1*h_2s2 - r2*(c11*gam(i, j, 3) + 2*c12*gam(i, j, 4) + c22*gam(i, j, 5));
            break;
        }
        case BC::None:
            break;
    }
}

void Membrane::n11BoundaryNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                           const double h_11s, const double z_1, const double z_2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            b(k) = (n11(i, j) - val/pow(h_11s, 2))/D;
            A.row(k).cols(    0,   nxy-1) = (H1111(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1112(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD1.row(k)
                                          + (H1112(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1122(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k).cols(  nxy, 2*nxy-1) = (H1111(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1112(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD1.row(k)
                                          + (H1112(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1122(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k).cols(2*nxy, 3*nxy-1) = H1111(i, j)*z_1*DD1.row(k) + H1112(i, j)*(z_1*DD2.row(k) + z_2*DD1.row(k)) + H1122(i, j)*z_2*DD2.row(k);
            A(k, k)     -= gam(i, j, 0)*(H1111(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1112(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 1)*(H1111(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1112(i, j)*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + H1122(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 2)*(H1112(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1122(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            A(k, k+nxy) -= gam(i, j, 3)*(H1111(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1112(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 4)*(H1111(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1112(i, j)*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + H1122(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 5)*(H1112(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1122(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            break;
        }
        case BC::Neumann:
            std::println("n11 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::Robin:
            std::println("n11 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n12BoundaryNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                    const double h_11s, const double h_12s, const double h_22s, const double z_1, const double z_2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            b(k+nxy) = (h_11s*h_12s*n11(i, j) + h_11s*h_22s*n12(i, j) - val)/D;
            double c11 = h_11s*h_12s*H1111(i, j) + h_11s*h_22s*H1112(i, j);
            double c12 = h_11s*h_12s*H1112(i, j) + h_11s*h_22s*H1221(i, j);
            double c22 = h_11s*h_12s*H1122(i, j) + h_11s*h_22s*H1222(i, j);
            A.row(k+nxy).cols(    0,   nxy-1) = (c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD1.row(k)
                                              + (c12*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k+nxy).cols(  nxy, 2*nxy-1) = (c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD1.row(k)
                                              + (c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k+nxy).cols(2*nxy, 3*nxy-1) = c11*z_1*DD1.row(k) + c12*(z_1*DD2.row(k) + z_2*DD1.row(k)) + c22*z_2*DD2.row(k);
            A(k+nxy, k)     -= gam(i, j, 0)*(c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 1)*(c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 2)*(c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            A(k+nxy, k+nxy) -= gam(i, j, 3)*(c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 4)*(c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 5)*(c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            break;
        }
        case BC::Neumann:
            std::println("n12 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::Robin:
            std::println("n12 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n21BoundaryNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                    const double h_11s, const double h_21s, const double h_22s, const double z_1, const double z_2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            b(k) = (h_11s*h_22s*n12(i, j) + h_21s*h_22s*n22(i, j) - val)/D;
            double c11 = h_11s*h_22s*H1112(i, j) + h_21s*h_22s*H1122(i, j);
            double c12 = h_11s*h_22s*H1221(i, j) + h_21s*h_22s*H1222(i, j);
            double c22 = h_11s*h_22s*H1222(i, j) + h_21s*h_22s*H2222(i, j);
            A.row(k).cols(    0,   nxy-1) = (c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD1.row(k)
                                          + (c12*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k).cols(  nxy, 2*nxy-1) = (c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD1.row(k)
                                          + (c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k).cols(2*nxy, 3*nxy-1) = c11*z_1*DD1.row(k) + c12*(z_1*DD2.row(k) + z_2*DD1.row(k)) + c22*z_2*DD2.row(k);
            A(k, k)     -= gam(i, j, 0)*(c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 1)*(c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 2)*(c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            A(k, k+nxy) -= gam(i, j, 3)*(c11*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + c12*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 4)*(c11*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c12*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + c22*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                         + gam(i, j, 5)*(c12*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + c22*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            break;
        }
        case BC::Neumann:
            std::println("n21 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
            break;
        case BC::Robin:
            std::println("n21 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

void Membrane::n22BoundaryNonlinear(const BC bc, const double val, const size_t i, const size_t j, arma::mat &A,
                                    const double h_22s, const double z_1, const double z_2)
{
    size_t k = i + j*nx;
    switch (bc)
    {
        case BC::Dirichlet:
        {
            b(k+nxy) = (n22(i, j) - val/pow(h_22s, 2))/D;
            A.row(k+nxy).cols(    0,   nxy-1) = (H1122(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD1.row(k)
                                              + (H1222(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H2222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k+nxy).cols(  nxy, 2*nxy-1) = (H1122(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1222(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD1.row(k)
                                              + (H1222(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H2222(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)))*DD2.row(k);
            A.row(k+nxy).cols(2*nxy, 3*nxy-1) = H1122(i, j)*z_1*DD1.row(k) + H1222(i, j)/2*(z_1*DD2.row(k) + z_2*DD1.row(k)) + H2222(i, j)*z_2*DD2.row(k);
            A(k+nxy, k)     -= gam(i, j, 0)*(H1122(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 1)*(H1122(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1222(i, j)*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + H2222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 2)*(H1222(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H2222(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            A(k+nxy, k+nxy) -= gam(i, j, 3)*(H1122(i, j)*(1 + e11(i, j)*v1__1(i, j) + e12(i, j)*v2__1(i, j)) + H1222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 4)*(H1122(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H1222(i, j)*(2 + e11(i, j)*v1__1(i, j) + e12(i, j)*(v2__1(i, j) + v1__2(i, j)) + e22(i, j)*v2__2(i, j)) + H2222(i, j)*(e11(i, j)*v1__2(i, j) + e12(i, j)*v2__2(i, j)))
                             + gam(i, j, 5)*(H1222(i, j)*(e12(i, j)*v1__1(i, j) + e22(i, j)*v2__1(i, j)) + H2222(i, j)*(1 + e12(i, j)*v1__2(i, j) + e22(i, j)*v2__2(i, j)));
            break;
        }
        case BC::Neumann:
            std::println("n22 Neumann not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::Robin:
            std::println("n22 Robin not implemented for planeStrain!");
            exit(EXIT_FAILURE);
        case BC::None:
            break;
    }
}

template void Membrane::boundary<int>(const Field, const Direction, const BC, const int);
template void Membrane::boundary<size_t>(const Field, const Direction, const BC, const size_t);
template void Membrane::boundary<double>(const Field, const Direction, const BC, const double);
template void Membrane::boundary<arma::vec>(const Field, const Direction, const BC, const arma::vec);
template void Membrane::boundary<int>(const Field, const Direction, const BC, const double, const double, const int);
template void Membrane::boundary<size_t>(const Field, const Direction, const BC, const double, const double, const size_t);
template void Membrane::boundary<double>(const Field, const Direction, const BC, const double, const double, const double);
template void Membrane::boundary<arma::vec>(const Field, const Direction, const BC, const double, const double, const arma::vec);
template void Membrane::boundary<int>(const Field, const Lagrange::CurveInterpolant*, const BC, const int);
template void Membrane::boundary<size_t>(const Field, const Lagrange::CurveInterpolant*, const BC, const size_t);
template void Membrane::boundary<double>(const Field, const Lagrange::CurveInterpolant*, const BC, const double);
template void Membrane::boundary<arma::vec>(const Field, const Lagrange::CurveInterpolant*, const BC, const arma::vec);
template void Membrane::boundary<int>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const int);
template void Membrane::boundary<size_t>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const size_t);
template void Membrane::boundary<double>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const double);
template void Membrane::boundary<arma::vec>(const Field, const Lagrange::CurveInterpolant*, const BC, const double, const double, const arma::vec);