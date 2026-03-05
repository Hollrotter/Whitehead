#include "VLM.hpp"

/**
 * @brief 
 * 
 */
void VLM::aerodynamicMatrix()
{
    // Needed for linear analysis only.
    // For the nonlinear analysis, the wake does make the
    // matrix dependent on alpha.
    nC = normalise(RA.cols(arma::uvec{2, 1})-RB.cols(arma::uvec{2, 1}), 2, 1);
    nC.col(1) *=-1;
    #pragma omp parallel for
    for (size_t k = 0; k < nxy; k++)
        A.col(k) = sum(horseshoe(RA.row(k), RB.row(k)).cols(1, 2)%nC, 1);
    if (sym == Symmetry::y)
    {
        RA.col(1) *=-1;
        RB.col(1) *=-1;
        #pragma omp parallel for
        for (size_t k = 0; k < nxy; k++)
            A.col(k) += sum(horseshoe(RB.row(k), RA.row(k)).cols(1, 2)%nC, 1);
        RA.col(1) *=-1;
        RB.col(1) *=-1;
    }
}