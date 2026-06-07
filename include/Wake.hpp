#pragma once
#include "Lagrange.hpp"

class Wake
{
    Lagrange::CurveInterpolant* chi;
public:
    Wake() = default;
    Wake(Lagrange::CurveInterpolant* _chi) : chi(_chi) {};
    friend class Wing;
    friend class Aerodynamics;
};