#pragma once
#include "misc.hpp"

class Interface
{
    size_t sourceDomain;
    size_t targetDomain;
    size_t sourceCurve;
    size_t targetCurve;
public:
    Interface(size_t sD, size_t tD, size_t sC, size_t tC) : sourceDomain(sD), targetDomain(tD), sourceCurve(sC), targetCurve(tC) {};
    friend class Structure;
};