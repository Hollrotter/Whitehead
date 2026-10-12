#include "MemDVM.hpp"

void MemDVM::linear()
{
    dvm->geometry();
    arma::mat TA = Lagrange::interpolationMatrix(string->x, dvm->xC);
    arma::mat TS = Lagrange::interpolationMatrix(dvm->xg, string->x);
    switch (coupling)
    {
        case Coupling::monolithic:
        {
            string->structuralMatrix();
            dvm->aerodynamicMatrix();
            size_t n_s = string->n;
            size_t n_a = dvm->nx;
            size_t n = n_s + n_a;
            arma::mat M = join_vert(join_horiz(     string->S, 2*dvm->qdyn*TS/string->sigma/repelem(arma::diff(dvm->x), 1, n_s).t()),
                                    join_horiz(-TA*string->D1, dvm->A));
            arma::mat RHS(n, dvm->con,  arma::fill::zeros);
            switch(string->front)
            {
                case BC::Dirichlet:
                    M.row(0).eye();
                    break;
                case BC::Neumann:
                    M.row(0) = string->D1.row(0);
                    break;
                case BC::Robin:
                    printf("Robin BC not implemented for String");
                    exit(EXIT_FAILURE);
                case BC::None:
                    printf("No BC provided at front");
                    exit(EXIT_FAILURE);
            }
            switch(string->back)
            {
                case BC::Dirichlet:
                    M.row(n_s-1).zeros();
                    M(n_s-1, n_s-1) = 1;
                    break;
                case BC::Neumann:
                    M.row(n_s-1) = string->D1.row(n_s-1);
                    break;
                case BC::Robin:
                    printf("Robin BC not implemented for String");
                    exit(EXIT_FAILURE);
                case BC::None:
                    printf("No BC provided at front");
                    exit(EXIT_FAILURE);
            }
            RHS.row(0).fill(string->frontBC);
            RHS.row(n_s-1).fill(string->backBC);
            arma::vec w = dvm->camber.diff(dvm->xC/dvm->c);
            RHS.rows(n_s, n-1) = repelem(w, 1, dvm->con) - repelem(dvm->alpha.t(), n_a, 1);
            RHS = solve(M, RHS);
            string->z   = RHS(arma::span(0, n_s-1), 0);
            arma::mat g = RHS.rows(n_s, n-1);
            dvm->postprocessing(g);
            break;
        }
        case Coupling::partitioned:
        {
            dvm->dvmSolve();
            string->solve_S();
            double L_old = 0, M_old = 0;
            for (size_t i = 1; i <= iter; i++)
            {
                std::println("Iteration {}/{}", i, iter);
                dvm->wE = TA*string->D1 * string->z;
                dvm->dvmEval();
                string->p = dvm->qdyn*TS * dvm->dcp.col(0);
                string->solve_b();
                double L_new = dvm->cL(0);
                double M_new = dvm->cM(0);
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