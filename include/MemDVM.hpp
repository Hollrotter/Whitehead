#pragma once
#include "String.hpp"
#include "DVM.hpp"
#include "Lagrange.hpp"

class MemDVM
{
    String* string; // Structure model
    DVM* dvm; // Aerodynamic model
    Coupling coupling = Coupling::monolithic; // Coupling approach
    size_t iter = 100; // Max. number of itertions
    double changeTarget = 1e-5; // Target change for cL and cM
public:
    MemDVM(String* s, DVM* d) : string(s), dvm(d) {};
    MemDVM(DVM* d, String* s) : string(s), dvm(d) {};
    void operator()(const Coupling c)
    {
        coupling = c;
    }
    // Linear analysis
    void linear();
};