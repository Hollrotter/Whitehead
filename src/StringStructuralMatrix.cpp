#include "String.hpp"

/**
 * @brief 
 * 
 */
void String::structuralMatrix()
{
    switch(analysis)
    {
        case Analysis::linear:
            S = D11;
            break;
        case Analysis::nonlinear:
            double s = sigma;
            arma::vec dzdx = D1*z;
            if (materialModel == Material::extensible)
            {
                double L = integrate(sqrt(1 + pow(dzdx, 2)));
                double delta = log(L/L0);
                s += Et*delta;
            }
            S = s*D11;
            b = s*D11*z + p%pow(1 + pow(dzdx, 2), 1.5);
            break;
    }
}