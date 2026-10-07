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
            membrane->structuralMatrix();

            arma::vec DX(nxy_v);
            for (size_t n = 0, k = 0; n < ny_v; n++)
                for (size_t m = 0; m < nx_v; m++, k++)
                    DX(k) = (vlm->x(m+1, n) - vlm->x(m, n) + vlm->x(m+1, n+1) - vlm->x(m, n+1))/2;
            vlm->aerodynamicMatrix();

            arma::mat M = join_vert(join_horiz(membrane->S, 2*vlm->qdyn*TS/repelem(DX, 1, nxy_m).t()),
                                    join_horiz(-TA * (repelem(vectorise(J_inv.slice(0)), 1, nxy_m) % membrane->DD1
                                                   +  repelem(vectorise(J_inv.slice(1)), 1, nxy_m) % membrane->DD2), vlm->A));

            arma::mat RHS(M.n_rows, vlm->con);
            arma::mat DDX = membrane->DD1;
            arma::mat DDY = membrane->DD2;
            for (size_t j = 0; j < ny_m; j++)
            {
                size_t k = j*nx_m;

                switch (membrane->z.southBC)
                {
                    case BC::Dirichlet:
                        M.row(k).zeros();
                        M(k, k) = 1;
                        RHS(k) = membrane->z.south(j);
                        break;
                    case BC::Neumann:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = DDX.row(k)*membrane->ec(0, j, 0) + DDY.row(k)*membrane->ec(0, j, 1);
                        RHS(k) = membrane->z.south(j);
                        break;
                    case BC::Robin:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = membrane->z.r2South*(DDX.row(k)*membrane->ec(0, j, 0) + DDY.row(k)*membrane->ec(0, j, 1));
                        M(k, k) += membrane->z.r1South;
                        RHS(k) = membrane->z.south(j);
                        break;
                    case BC::None:
                        break;
                    default:
                        std::println("A boundary condition was chosen, that is not implemented for membranes!");
                        exit(EXIT_FAILURE);
                }
                size_t i = nx_m-1;
                k = j*nx_m + i;
                switch (membrane->z.northBC)
                {
                    case BC::Dirichlet:
                        M.row(k).zeros();
                        M(k, k) = 1;
                        RHS(k) = membrane->z.north(j);
                        break;
                    case BC::Neumann:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = DDX.row(k)*membrane->ec(i, j, 0) + DDY.row(k)*membrane->ec(i, j, 1);
                        RHS(k) = membrane->z.north(j);
                        break;
                    case BC::Robin:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = membrane->z.r2North*(DDX.row(k)*membrane->ec(i, j, 0) + DDY.row(k)*membrane->ec(i, j, 1));
                        M(k, k) += membrane->z.r1North;
                        RHS(k) = membrane->z.north(j);
                        break;
                    case BC::None:
                        break;
                    default:
                        std::println("A boundary condition was chosen, that is not implemented for membranes!");
                        exit(EXIT_FAILURE);
                }
            }
            for (size_t i = 0; i < nx_m; i++)
            {
                switch (membrane->z.westBC)
                {
                    case BC::Dirichlet:
                        M.row(i).zeros();
                        M(i, i) = 1;
                        RHS(i) = membrane->z.west(i);
                        break;
                    case BC::Neumann:
                        M.row(i).zeros();
                        M.row(i).cols(0, nxy_m-1) = DDY.row(i)*membrane->ec(i, 0, 2) + DDX.row(i)*membrane->ec(i, 0, 1);
                        RHS(i) = membrane->z.west(i);
                        break;
                    case BC::Robin:
                        M.row(i).zeros();
                        M.row(i).cols(0, nxy_m-1) = membrane->z.r2West*(DDY.row(i)*membrane->ec(i, 0, 2) + DDX.row(i)*membrane->ec(i, 0, 1));
                        M(i, i) += membrane->z.r1West;
                        RHS(i) = membrane->z.west(i);
                        break;
                    case BC::None:
                        break;
                    default:
                        std::println("A boundary condition was chosen, that is not implemented for membranes!");
                        exit(EXIT_FAILURE);
                }
                size_t j = ny_m-1;
                size_t k = i + j*nx_m;
                switch (membrane->z.eastBC)
                {
                    case BC::Dirichlet:
                        M.row(k).zeros();
                        M(k, k) = 1;
                        RHS(k) = membrane->z.east(i);
                        break;
                    case BC::Neumann:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = DDY.row(i)*membrane->ec(i, j, 2) + DDX.row(i)*membrane->ec(i, j, 1);
                        RHS(k) = membrane->z.east(i);
                        break;
                    case BC::Robin:
                        M.row(k).zeros();
                        M.row(k).cols(0, nxy_m-1) = membrane->z.r2East*(DDY.row(i)*membrane->ec(i, j, 2) + DDX.row(i)*membrane->ec(i, j, 1));
                        M(k, k) += membrane->z.r1East;
                        RHS(k) = membrane->z.east(i);
                        break;
                    case BC::None:
                        break;
                    default:
                        std::println("A boundary condition was chosen, that is not implemented for membranes!");
                        exit(EXIT_FAILURE);
                }
            }
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