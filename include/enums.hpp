#pragma once

// Analysis type (linear/nonlinear)
enum class Analysis
{
    linear, // Linear analysis
    semilinear, // Semilinear analysis
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

// Describes the location of the wing boundary (north, south, east, west)
enum class Direction
{
    S, // Southern boundary (Normally the leading edge)
    E, // Eastern boundary (Normally the tip of the wing)
    N, // Northern boundary (Normally the trailing edge)
    W // Western boundary (Normally the base of the wing)
};

// For output or boundary conditions, one field must be selected
enum class Field
{
    z, // Deformation in z-direction
    v1, // Displacement in 1-direction (physical component)
    v2, // Displacement in 2-direction (physical component)
    n11, // Stress component 11 (physical component)
    n12, // Stress component 12 (physical component)
    n22  // Stress component 22 (physical component)
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