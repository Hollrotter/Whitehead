#pragma once
#include "Membrane.hpp"
#include "VLM.hpp"

class MemVLM
{
    Membrane* membrane; // Structure model
    VLM* vlm; // Aerodynamic model
    Coupling coupling = Coupling::monolithic; // Coupling approach
    size_t iter = 100; // Max. number of iterations for partitioned solution
    double changeTarget = 1e-5; // Target change for cL and cM in partitioned solution
public:
    MemVLM(Membrane* m, VLM* v) : membrane(m), vlm(v) {};
    MemVLM(VLM* v, Membrane* m) : membrane(m), vlm(v) {};
    // Linear analysis
    void linear();
};