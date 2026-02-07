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
        case BC::None:
            std::println("No BC provided at front!");
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
        case BC::None:
            std::println("No BC provided at front!");
            exit(EXIT_FAILURE);
    }
    arma::lu(L, U, P, S);
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
    z = solve(trimatu(U), solve(trimatl(L), P*b));
}