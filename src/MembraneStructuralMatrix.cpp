#include "Membrane.hpp"

/**
 * @brief 
 * 
 */
void Membrane::structuralMatrix()
{
    if (analysis == Analysis::linear || analysis == Analysis::semilinear)
        S = d2dx2(n11) + d2dxdy(2*n12) + d2dy2(n22)
        - ddx(gam.slice(0)%n11 + 2*gam.slice(1)%n12 + gam.slice(2)%n22)
        - ddy(gam.slice(3)%n11 + 2*gam.slice(4)%n12 + gam.slice(5)%n22);
    else
        std::println("structuralMatrix is currently not needed for nonlinear!");
}