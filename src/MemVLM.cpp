#include "MemVLM.hpp"

void MemVLM::linear()
{
    vlm->geometry();
    size_t nx_v = vlm->nx;
    size_t ny_v = vlm->ny;
    size_t nx_m = membrane->nx;
    size_t ny_m = membrane->ny;

    std::println("Creating interpolation matrix for displacement");
    arma::mat  x_m_hat = Chebyshev::DiscreteChebyshevTransform(membrane->x1, membrane->x2, membrane->x);
    arma::mat  y_m_hat = Chebyshev::DiscreteChebyshevTransform(membrane->x1, membrane->x2, membrane->y);
    auto [x1_m_bar, x2_m_bar] = Chebyshev::DerivativeCoefficients(x_m_hat);
    auto [y1_m_bar, y2_m_bar] = Chebyshev::DerivativeCoefficients(y_m_hat);

    arma::mat xC = repelem(Chebyshev::gaussLobatto(nx_v),     1, ny_v);
    arma::mat yC = repelem(Chebyshev::gaussLobatto(ny_v).t(), nx_v, 1);
    
    #pragma omp parallel for
    for (size_t i = 0; i < nx_v; i++)
        for (size_t j = 0; j < ny_v; j++)
            for (size_t iter = 0; iter < 50; iter++)
            {
                arma::mat J(2, 2, arma::fill::zeros);
                arma::vec b = {vlm->rC(i, j, 0), vlm->rC(i, j, 1)};
                for (size_t p = 0; p < nx_m; p++)
                {
                    double Tx = boost::math::chebyshev_t(p,-xC(i, j));
                    for (size_t q = 0; q < ny_m; q++)
                    {
                        double Txy = Tx*boost::math::chebyshev_t(q,-yC(i, j));
                        J += Txy*arma::mat{{x1_m_bar(p, q), x2_m_bar(p, q)},
                                           {y1_m_bar(p, q), y2_m_bar(p, q)}};
                        b -= Txy*arma::vec{x_m_hat(p, q), y_m_hat(p, q)};
                    }
                }
                b = solve(J, b);
                xC(i, j) -= 0.2*b(0);
                yC(i, j) -= 0.2*b(1);
                if (norm(b) < 1e-6)
                    break;
            }
    arma::mat TA = Lagrange::interpolationMatrix(membrane->x1, membrane->x2, xC, yC);

    std::println("Creating interpolation matrix for pressure");
    arma::vec xG = Chebyshev::gaussLobatto(nx_v);
    arma::vec yG = Chebyshev::gaussLobatto(ny_v);
    arma::mat  x_w_hat = Chebyshev::DiscreteChebyshevTransform(xG, yG, vlm->rG.slice(0));
    arma::vec  y_w_hat = Chebyshev::DiscreteChebyshevTransform(yG, vlm->rG.slice(1).row(0).t());
    auto [x1_w_bar, x2_w_bar] = Chebyshev::DerivativeCoefficients(x_w_hat);
    arma::vec  y_w_bar = Chebyshev::DerivativeCoefficients(y_w_hat);

    arma::mat X_m = repelem(membrane->x1,     1, ny_m);
    arma::mat Y_m = repelem(membrane->x2.t(), nx_m, 1);

    #pragma omp parallel for
    for (size_t i = 0; i < nx_m; i++)
        for (size_t j = 0; j < ny_m; j++)
            for (size_t iter = 0; iter < 100; iter++)
            {
                arma::mat J(2, 2, arma::fill::zeros);
                arma::vec b = {membrane->x(i, j), membrane->y(i, j)};
                for (size_t q = 0; q < ny_v; q++)
                {
                    double Ty = boost::math::chebyshev_t(q,-Y_m(i, j));
                    J(1, 1) +=  Ty*y_w_bar(q);
                    b(1)    -=  Ty*y_w_hat(q);
                    for (size_t p = 0; p < nx_v; p++)
                    {
                        double Txy = boost::math::chebyshev_t(p,-X_m(i, j))*Ty;
                        J(0, 0) += Txy*x1_w_bar(p, q);
                        J(0, 1) += Txy*x2_w_bar(p, q);
                        b(0)    -= Txy*x_w_hat(p, q);
                    }
                }
                b = solve(trimatu(J), b);
                X_m(i, j) -= 0.2*b(0);
                Y_m(i, j) -= 0.2*b(1);
                if (norm(b) < 1e-10)
                    break;
            }
    arma::mat TS = Lagrange::interpolationMatrix(xG, yG, X_m, Y_m);

    arma::mat det = membrane->J(0, 0)%membrane->J(1, 1) - membrane->J(0, 1)%membrane->J(1, 0);
    arma::cube J_inv = join_slices(membrane->J(1, 1)/det, -membrane->J(0, 1)/det);

    switch (coupling)
    {
        case Coupling::monolithic:
        {
            size_t nxy_m = membrane->nxy;
            size_t nxy_v = vlm->nxy;
            size_t nxy = nxy_m + nxy_v;
            membrane->structuralMatrix();
            for (size_t j = 1; j < ny_m-1; j++)
            {
                // BC west
                membrane->zBoundary(membrane->z.westBC, membrane->z.west(j),      0, j, membrane->z.r1West, membrane->z.r2West, membrane->h_1s1_west(j), membrane->h_1s2_west(j));
                // BC east
                membrane->zBoundary(membrane->z.eastBC, membrane->z.east(j), nx_m-1, j, membrane->z.r1East, membrane->z.r2East, membrane->h_1s1_east(j), membrane->h_1s2_east(j));
            }

            for (size_t i = 1; i < nx_m-1; i++)
            {
                // BC south
                membrane->zBoundary(membrane->z.southBC, membrane->z.south(i), i,      0, membrane->z.r1South, membrane->z.r2South, membrane->h_2s1_south(i), membrane->h_2s2_south(i));
                // BC north
                membrane->zBoundary(membrane->z.northBC, membrane->z.north(i), i, ny_m-1, membrane->z.r1North, membrane->z.r2North, membrane->h_2s1_north(i), membrane->h_2s2_north(i));
            }

            // BC south-west corner (i = 0, j = 0)
            membrane->zBoundary(membrane->z.southBC, membrane->z.south(0), 0, 0, membrane->z.r1South, membrane->z.r2South, membrane->h_2s1_south(0), membrane->h_2s2_south(0));

            // BC north-west corner (i = 0, j = ny-1)
            size_t j = ny_m-1;
            membrane->zBoundary(membrane->z.westBC, membrane->z.west(j),   0, j, membrane->z.r1West,  membrane->z.r2West,  membrane->h_1s1_west(j),  membrane->h_1s2_west(j));

            // BC south-east corner (i = nx-1, j = 0)
            size_t i = nx_m-1;
            membrane->zBoundary(membrane->z.eastBC, membrane->z.east(0),   i, 0, membrane->z.r1East,  membrane->z.r2East,  membrane->h_1s1_east(0),  membrane->h_1s2_east(0));

            // BC north-east corner (i = nx-1, j = ny-1)
            i = nx_m-1;
            j = ny_m-1;
            membrane->zBoundary(membrane->z.northBC, membrane->z.north(i), i, j, membrane->z.r1North, membrane->z.r2North, membrane->h_2s1_north(i), membrane->h_2s2_north(i));

            arma::vec DX(nxy_v);
            for (size_t n = 0, k = 0; n < ny_v; n++)
                for (size_t m = 0; m < nx_v; m++, k++)
                    DX(k) = (vlm->x(m+1, n) - vlm->x(m, n) + vlm->x(m+1, n+1) - vlm->x(m, n+1))/2;
            vlm->aerodynamicMatrix();

            arma::mat M = join_vert(join_horiz(membrane->S, 2*vlm->qdyn*TS/repelem(DX, 1, nxy_m).t()),
                                    join_horiz(-TA * (repelem(vectorise(J_inv.slice(0)), 1, nxy_m) % membrane->DD1
                                                   +  repelem(vectorise(J_inv.slice(1)), 1, nxy_m) % membrane->DD2), vlm->A));

            for (size_t j = 0; j < ny_m; j++)
            {
                size_t k = j*nx_m;
                M(k, arma::span(nxy_m, nxy-1)).zeros();

                size_t i = nx_m-1;
                k = j*nx_m + i;
                M(k, arma::span(nxy_m, nxy-1)).zeros();
            }

            for (size_t i = 0; i < nx_m; i++)
            {
                M(i, arma::span(nxy_m, nxy-1)).zeros();

                size_t j = ny_m-1;
                size_t k = i + j*nx_m;
                M(k, arma::span(nxy_m, nxy-1)).zeros();
            }

            arma::mat RHS(M.n_rows, vlm->con);
            for (size_t c = 0; c < vlm->con; c++)
                RHS(arma::span(0, nxy_m-1), c) = membrane->b(arma::span(0, nxy_m-1));

            arma::vec r = repelem((vlm->ar.head(ny_v)+vlm->ar.tail(ny_v))/2, nx_v, 1);
            arma::vec w = arma::zeros(nxy_v);
            #pragma omp parallel for
            for (size_t n = 0; n < ny_v; n++)
                w.subvec(n*nx_v, (n+1)*nx_v-1)
                    = vlm->c.diff((2*vlm->RC(arma::span(n*nx_v,(n+1)*nx_v-1),0)-vlm->x(0,n)-vlm->x(0,n+1))
                    /(vlm->x(nx_v,n)-vlm->x(0,n)+vlm->x(nx_v,n+1)-vlm->x(0,n+1)));
            RHS.rows(nxy_m, M.n_rows-1) = repelem(w - r, 1, vlm->con) - vlm->nC.col(1)*vlm->alpha.t();
            RHS = solve(M, RHS);
            membrane->z = arma::reshape(RHS(arma::span(0, nxy_m-1), 0), nx_m, ny_m);
            arma::mat g = RHS.rows(nxy_m, M.n_rows-1);
            vlm->postprocessing(g);
            break;
        }
        case Coupling::partitioned:
        {
            std::println("Solving VLM");
            vlm->vlmSolve();
            std::println("Solving membrane");
            membrane->solve_S();

            double L_old = 0;
            double M_old = 0;
            for (size_t n = 1; n <= iter; n++)
            {
                std::println("Iteration {}/{}", n, iter);

                vlm->wE = TA*vectorise((membrane->D1*membrane->z)%J_inv.slice(0) + (membrane->z*membrane->D2.t())%J_inv.slice(1));
                vlm->vlmEval();

                membrane->p = arma::reshape(vlm->qdyn*TS*vectorise(vlm->dcp), nx_m, ny_m);
                membrane->solve_b();

                double L_new = vlm->lift(0);
                double M_new = vlm->moment(0);
                double changeL = fabs(1-L_old/L_new);
                double changeM = fabs(1-M_old/M_new);
                std::println("Relative change of Lift:   {0:4.2e}",   changeL);
                std::println("Relative change of Moment: {0:4.2e}\n", changeM);
                if (changeL < changeTarget && changeM < changeTarget)
                    break;
                L_old = L_new;
                M_old = M_new;
            }
            break;
        }
    }
}