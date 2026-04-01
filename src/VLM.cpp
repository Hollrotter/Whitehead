#include "VLM.hpp"

/**
 * @brief 
 * 
 * @param _qdyn Dynamic pressure of the inflow.
 */
void VLM::dynamicPressure(double _qdyn)
{
    qdyn = _qdyn;
    try
    {
        if (qdyn <= 0)
            throw std::runtime_error("Dynamic pressure must be positive!");
    }
    catch(const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief 
 * 
 * @param _alpha Pitch given as vector for multiple configurations.
 */
void VLM::pitch(arma::vec _alpha)
{
    alpha = arma::datum::pi/180*_alpha;
    con = alpha.size();
    lift.zeros(con);
    moment.zeros(con);
    dcp.zeros(nx, ny, con);
}

void VLM::vlm()
{
    switch(analysis)
    {
        case Analysis::linear:
            vlmSolve();
            vlmEval();
            break;
        case Analysis::nonlinear:
            vlmNonlinear();
            break;
    }
}

/**
 * @brief 
 * 
 * @param filename Path for the file data will be written to.
 */
void VLM::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t m = 0; m < nx; m++, file << '\n')
        for (size_t n = 0; n < ny; n++, file << '\n')
        {
            file << rG(m, n, 0) << ' ' << rG(m, n, 1);
            for (size_t i = 0; i < con; i++)
                file << ' ' << dcp(m, n, i);
        }
    file.close();
}

/**
 * @brief 
 * 
 */
void VLM::geometry()
{
    arma::mat x00 = x.submat(0, 0, nx-1, ny-1);
    arma::mat x10 = x.submat(1, 0, nx,   ny-1);
    arma::mat x01 = x.submat(0, 1, nx-1, ny);
    arma::mat x11 = x.submat(1, 1, nx,   ny);
    arma::mat z00 = z.submat(0, 0, nx-1, ny-1);
    arma::mat z10 = z.submat(1, 0, nx,   ny-1);
    arma::mat z01 = z.submat(0, 1, nx-1, ny);
    arma::mat z11 = z.submat(1, 1, nx,   ny);
    RA.col(0) = vectorise((3*x00 + x10)/4);
    RA.col(1) = repelem(y.head(ny), nx, 1);
    RA.col(2) = vectorise((3*z00 + z10)/4);
    RB.col(0) = vectorise((3*x01 + x11)/4);
    RB.col(1) = repelem(y.tail(ny), nx, 1);
    RB.col(2) = vectorise((3*z01 + z11)/4);
    RC.col(0) = vectorise((3*(x10 + x11) + x00 + x01)/8);
    RC.cols(1, 2) = (RA.cols(1, 2)+RB.cols(1, 2))/2;

    #pragma omp parallel for
    for (size_t d = 0; d < 3; d++)
    {
        rC.slice(d) = reshape(RC.col(d), nx, ny);
        rG.slice(d) = reshape((RA.col(d)+RB.col(d))/2, nx, ny);
    }
}

void VLM::vlmSolve()
{
    z = repelem(z0.t(), nx+1, 1);
    geometry();
    aerodynamicMatrix();
    arma::lu(L, U, P, A);
}

void VLM::vlmEval()
{
    arma::vec r = repelem((ar.head(ny)+ar.tail(ny))/2, nx, 1);
    arma::vec w = arma::zeros(nx*ny);
    #pragma omp parallel for
    for (size_t n = 0; n < ny; n++)
        w.subvec(n*nx,(n+1)*nx-1) += c.diff((2*RC(arma::span(n*nx,(n+1)*nx-1),0)-x(0,n)-x(0,n+1))/(x(nx, n)-x(0, n)+x(nx, n+1)-x(0, n+1)));
    arma::mat g = solve(trimatu(U), solve(trimatl(L), P*(repelem(w - r, 1, con) - nC.col(1)*alpha.t())));
    postprocessing(g);
}

void VLM::vlmNonlinear()
{
    #pragma omp parallel for
    for (size_t n = 0; n < ny+1; n++)
    {
        double chord = x(nx, n) - x(0, n);
        if (almostEqual(chord, 0))
            z.col(n).fill(z0(n));
        else
            for (size_t m = 0; m < nx+1; m++)
                z(m, n) = z0(n) + chord*c((x(m, n)-x(0, n))/chord);
    }
    geometry();

    // Normal calculated by modelling panel with bilinear function and analytically calculating the gradient
    nC.set_size(nxy, 3);
    #pragma omp parallel for
    for (size_t n = 0; n < ny; n++)
        for (size_t m = 0; m < nx; m++)
        {
            size_t k = m + n*nx;
            double c1 = (z(m, n) - z(m+1, n))/(x(m, n) - x(m+1, n));
            double c2 = ((x(m, n) - x(m+1, n+1))/(x(m, n+1) - x(m+1, n+1))*z(m, n+1) - z(m, n) + (x(m, n+1) - x(m, n))/(x(m, n+1) - x(m+1, n+1))*z(m+1, n+1))/(y(n+1) - y(n));
            double c3 = ((z(m+1, n+1) - z(m, n+1))/(x(m+1, n+1) - x(m, n+1)) - (z(m+1, n) - z(m, n))/(x(m+1, n) - x(m, n)))/(y(n+1) - y(n));
            nC.row(k) = normalise(arma::rowvec{c3*(y(n)-RC(k, 1))-c1, c3*(x(m, n)-RC(k, 0))-c2, 1});
        }
    #pragma omp parallel for
    for (size_t n = 0; n < ny; n++)
        for (size_t m = 0; m < nx-1; m++)
        {
            size_t k = m + n*nx;
            arma::rowvec::fixed<3> r3 = {x(m+1, n+1), y(n+1), z(m+1, n+1)};
            arma::rowvec::fixed<3> r6 = {x(m+1,   n), y(n),   z(m+1,   n)};
            A.col(k) = sum(ring({RA.row(k), RB.row(k), r3, RB.row(k+1), RA.row(k+1), r6})%nC, 1);
        }
    if (sym == Symmetry::y)
    {
        RA.col(1) *=-1;
        RB.col(1) *=-1;
        #pragma omp parallel for
        for (size_t n = 0; n < ny; n++)
            for (size_t m = 0; m < nx-1; m++)
            {
                size_t k = m + n*nx;
                arma::rowvec::fixed<3> r3 = {x(m+1, n+1),-y(n+1), z(m+1, n+1)};
                arma::rowvec::fixed<3> r6 = {x(m+1,   n),-y(n),   z(m+1,   n)};
                A.col(k) += sum(ring({RB.row(k), RA.row(k), r6, RA.row(k+1), RB.row(k+1), r3})%nC, 1);
            }
        RA.col(1) *=-1;
        RB.col(1) *=-1;
    }
    arma::mat g(nxy, con, arma::fill::zeros);
    for (size_t i = 0; i < alpha.size(); i++)
    {
        #pragma omp parallel for
        for (size_t n = 0; n < ny; n++)
        {
            size_t k = nx*(n+1) - 1;
            arma::rowvec::fixed<3> r3 = (4*arma::rowvec::fixed<3>({x(nx, n+1), y(n+1), z(nx, n+1)}) - RB.row(k))/3;
            arma::rowvec::fixed<3> r6 = (4*arma::rowvec::fixed<3>({x(nx,   n), y(n),   z(nx,   n)}) - RA.row(k))/3;
            A.col(k) = sum((line(r6, RA.row(k)) + line(RA.row(k), RB.row(k)) + line(RB.row(k), r3) + wake(i, r6, r3))%nC, 1);
        }
        if (sym == Symmetry::y)
        {
            RA.col(1) *=-1;
            RB.col(1) *=-1;
            #pragma omp parallel for
            for (size_t n = 0; n < ny; n++)
            {
                size_t k = nx*(n+1) - 1;
                arma::rowvec::fixed<3> r3 = (4*arma::rowvec::fixed<3>({x(nx, n+1),-y(n+1), z(nx, n+1)}) - RB.row(k))/3;
                arma::rowvec::fixed<3> r6 = (4*arma::rowvec::fixed<3>({x(nx,   n),-y(n),   z(nx,   n)}) - RA.row(k))/3;
                A.col(k) += sum((line(r3, RB.row(k)) + line(RB.row(k), RA.row(k)) + line(RA.row(k), r6) + wake(i, r3, r6))%nC, 1);
            }
            RA.col(1) *=-1;
            RB.col(1) *=-1;
        }
        g.col(i) = solve(A, -nC*arma::vec::fixed<3>{cos(alpha(i)), 0, sin(alpha(i))});
    }
    postprocessing(g);
}

/**
 * @brief 
 * 
 * @param g Circulation/inflow velocity calculated by VLM.
 */
void VLM::postprocessing(arma::mat &g)
{
    area = 0;
    lift.zeros();
    moment.zeros();
            switch(analysis)
            {
                case Analysis::linear:
                    for (size_t n = 0, k = 0; n < ny; n++)
                    {
                        double dy = y(n+1) - y(n);
                        for (size_t m = 0; m < nx; m++, k++)
                        {
                            double dx = (x(m+1, n) - x(m, n) + x(m+1, n+1) - x(m, n+1))/2;
                            area += dx*dy;
                            dcp.tube(m, n) = 2/dx*g.row(k);
                            lift   += 2*dy*g.row(k).t();
                            moment -= 2*dy*g.row(k).t()*rG(m, n, 0);
                        }
                    }
                    break;
                case Analysis::nonlinear:
                    for (size_t n = 0, k = 0; n < ny; n++)
                    {
                        double dy = y(n+1) - y(n);
                        for (size_t m = 0; m < nx; m++, k++)
                        {
                            // double g_k = g(k, 0);
                            arma::vec::fixed<3> v1 = {x(m+1,   n) - x(m, n+1), y(n) - y(n+1), z(m+1,   n) - z(m, n+1)};
                            arma::vec::fixed<3> v2 = {x(m+1, n+1) - x(m,   n), y(n+1) - y(n), z(m+1, n+1) - z(m,   n)};
                            arma::vec::fixed<3> n_k = cross(v1, v2)/2;
                            double dA = norm(n_k)/2;
                            area += dA;
                            double dL = m > 0 ? (g(k, 0) - g(k-1, 0))*dy : g(k, 0)*dy;
                            lift   += 2*dL;
                            dcp(m, n, 0) = dL/dA;
                            // moment -= 2*Fz*rG(m, n, 0); // Evtl. nicht richtig!
                        }
                    }
                    break;
            }
            lift   *= qdyn;
            moment *= qdyn;
}