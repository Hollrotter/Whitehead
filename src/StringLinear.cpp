#include "String.hpp"

/**
 * @brief 
 * 
 */
void String::linear()
{
    solve_S();
    solve_b();
}

/**
 * @brief 
 * 
 */
void String::solve_S()
{
    analysis = Analysis::linear;
    structuralMatrix();
    switch(front)
    {
        case BC::Dirichlet:
            S.row(0).eye();
            break;
        case BC::Neumann:
            S.row(0) = D1.row(0);
            break;
        case BC::Robin:
            std::println("Robin BC not implemented for String!");
            exit(EXIT_FAILURE);
        default:
            std::println("No suitable BC provided at front!");
            exit(EXIT_FAILURE);
    }
    switch(back)
    {
        case BC::Dirichlet:
            S.row(n-1).zeros();
            S(n-1, n-1) = 1;
            break;
        case BC::Neumann:
            S.row(n-1) = D1.row(n-1);
            break;
        case BC::Robin:
            std::println("Robin BC not implemented for String!");
            exit(EXIT_FAILURE);
        default:
            std::println("No suitable BC provided at back!");
            exit(EXIT_FAILURE);
    }
    arma::mat P;
    arma::lu(L, U, P, S);
    for (size_t k = 0; k < n; k++)
        perm_idx(k) = P.row(k).index_max();
}

/**
 * @brief 
 * 
 */
void String::solve_b()
{
    b =-p/sigma;
    b(0)   = frontBC;
    b(n-1) = backBC;
    arma::vec b_perm = b.elem(perm_idx);
    arma::vec z_bar = solve(trimatl(L), b_perm, arma::solve_opts::fast);
    z               = solve(trimatu(U), z_bar,  arma::solve_opts::fast);
}