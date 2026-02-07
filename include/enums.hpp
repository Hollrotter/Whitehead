#pragma once

// Analysis type (linear/nonlinear)
enum class Analysis
{
    linear, // Linear analysis
    nonlinear // Nonlinear analysis
};

// Boundary type
enum class BC
{
    Dirichlet, // Boundary condition with predescribed value
    Neumann, // Boundary condition with predescribed derivative/gradient
    Robin, // Mixed boundary condition
    None // No boundary condition
};

// First or second derivative
enum class Derivative
{
    first, // First derivative
    second // Second derivative
};

// Boundary location of the camber (leading- or trailing edge)
enum class Location
{
    Front, // Front boundary (Normally the leading edge of the camber)
    Back // Back boundary (Normally the trailing edge of the camber)
};

// Material model (inextensible or extensible)
enum class Material
{
    inextensible, // No change of String length
    extensible // With change of String length
};