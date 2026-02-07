#include "String.hpp"

/**
 * @brief 
 * 
 */
void String::nonlinear()
{
    analysis = Analysis::nonlinear;
    for (size_t q = 0; q < iter; q++)
    {
        structuralMatrix();
        arma::vec dzdx = D1*z;
        arma::mat J = S + 3*((p%sqrt(1 + pow(dzdx, 2)))*dzdx.t())%D1;
        
        switch(front)
        {
            case BC::Dirichlet:
                J.row(0).eye();
                b(0) = z(0) - frontBC;
                break;
            case BC::Neumann:
                J.row(0) = D1.row(0);
                b(0) = dot(D1.row(0), z) - frontBC;
                break;
            case BC::Robin:
                std::println("Robin BC not implemented for String");
                exit(EXIT_FAILURE);
            case BC::None:
                std::println("No BC provided at front");
                exit(EXIT_FAILURE);
        }
        switch(back)
        {
            case BC::Dirichlet:
                J.row(n-1).zeros();
                J(n-1, n-1) = 1;
                b(n-1) = z(n-1) - backBC;
                break;
            case BC::Neumann:
                J.row(n-1) = D1.row(n-1);
                b(n-1) = backBC - dot(D1.row(n-1), z);
                break;
            case BC::Robin:
                std::println("Robin BC not implemented for String");
                exit(EXIT_FAILURE);
            case BC::None:
                std::println("No BC provided at front");
                exit(EXIT_FAILURE);
        }
            
        z -= omega*solve(J, b);
        double resMAX = max(arma::abs(b));
        double resRMS = sqrt(sum(b%b)/n);
        std::println("Iteration {}/{}", q+1, iter);
        std::println("MAX Residual: {}",   resMAX);
        std::println("RMS Residual: {}\n", resRMS);
        if (resRMS < residualTarget)
            break;
    }
}