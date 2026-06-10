#pragma once
#include "Lagrange.hpp"

class Wake
{
    Lagrange::CurveInterpolant* chi;
public:
    explicit Wake(Lagrange::CurveInterpolant* _chi) : chi(_chi) {};
    friend class Wing;
    friend class Aerodynamics;
};