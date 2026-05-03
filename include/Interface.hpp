#pragma once
#include "misc.hpp"

class Interface
{
    size_t sourceDomain;
    size_t targetDomain;
    size_t sourceCurve;
    size_t targetCurve;
    double lambdaSource = 0;
    double lambdaTarget = 0;
    arma::vec c1Source;
    arma::vec c1Target;
    arma::vec c2Source;
    arma::vec c2Target;
    arma::vec c11Source;
    arma::vec c11Target;
    arma::vec c12Source;
    arma::vec c12Target;
    arma::vec c22Source;
    arma::vec c22Target;
public:
    Interface(size_t sD, size_t tD, size_t sC, size_t tC) : sourceDomain(sD), targetDomain(tD), sourceCurve(sC), targetCurve(tC) {};
    friend class Structure;
    friend class Aerodynamics;
};