#include "VLM.hpp"

/**
 * @brief 
 * 
 * @param rA 
 * @param rB 
 * @return arma::mat 
 */
arma::mat VLM::line(arma::rowvec rA, arma::rowvec rB)
{
    arma::mat v(nxy, 3);
    #pragma omp parallel for
    for (size_t k = 0; k < nxy; k++)
    {
        arma::rowvec::fixed<3> r1 = RC.row(k) - rA;
        arma::rowvec::fixed<3> r2 = RC.row(k) - rB;
        double r1_n = norm(r1);
        double r2_n = norm(r2);
        if (r1_n < arma::datum::eps || r2_n < arma::datum::eps || r1_n*r2_n + dot(r1, r2) < arma::datum::eps)
            v.row(k).zeros();
        else
            v.row(k) = (r1_n + r2_n)*cross(r1, r2)/r1_n/r2_n/(r1_n*r2_n + dot(r1, r2));
    }
    return v/(4*arma::datum::pi);
}

/**
 * @brief 
 * 
 * @param rA 
 * @param rB 
 * @return arma::mat 
 */
arma::mat VLM::horseshoe(arma::rowvec rA, arma::rowvec rB)
{
    arma::mat v(nxy, 3);
    const arma::rowvec::fixed<3> ex = {1, 0, 0};
    arma::rowvec::fixed<3> rAB = rB - rA;
    #pragma omp parallel for
    for (size_t k = 0; k < nxy; k++)
    {
        arma::rowvec::fixed<3> rAP = RC.row(k) - rA;
        arma::rowvec::fixed<3> rBP = RC.row(k) - rB;
        arma::rowvec::fixed<3>  cA = cross( ex, rAP);
        arma::rowvec::fixed<3>  cP = cross(rAP, rBP);
        arma::rowvec::fixed<3>  cB = cross( ex, rBP);
        v.row(k) =-cA/dot(cA, cA)*(1 + dot(ex, normalise(rAP)))
                 + cP/dot(cP, cP)*(dot(rAB, normalise(rAP)) - dot(rAB, normalise(rBP)))
                 + cB/dot(cB, cB)*(1 + dot(ex, normalise(rBP)));
    }
    return v/(4*arma::datum::pi);
}

/**
 * @brief 
 * 
 * @param r 
 * @return arma::mat 
 */
arma::mat VLM::ring(arma::field<arma::rowvec> r)
{
    arma::mat v = line(r.back(), r(0));
    for (size_t i = 0; i < r.size()-1; i++)
        v += line(r(i), r(i+1));
    return v;
}

/**
 * @brief 
 * 
 * @param i 
 * @param rA 
 * @param rB 
 * @return arma::mat 
 */
arma::mat VLM::wake(size_t i, arma::rowvec rA, arma::rowvec rB)
{
    arma::mat v(nxy, 3);
    arma::rowvec::fixed<3> e = {cos(alpha(i)), 0, sin(alpha(i))};
    #pragma omp parallel for
    for (size_t k = 0; k < nxy; k++)
    {
        arma::rowvec::fixed<3> rAP = RC.row(k) - rA;
        arma::rowvec::fixed<3> rBP = RC.row(k) - rB;
        arma::rowvec::fixed<3>  cA = cross(e, rAP);
        arma::rowvec::fixed<3>  cB = cross(e, rBP);
        v.row(k) = cB/dot(cB, cB)*(1 + dot(e, normalise(rBP)))
                 - cA/dot(cA, cA)*(1 + dot(e, normalise(rAP)));
    }
    return v/(4*arma::datum::pi);
}