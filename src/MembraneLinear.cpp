#include "Membrane.hpp"

void Membrane::linear()
{
    analysis = Analysis::linear;
    solve_S();
    solve_b();
}

void Membrane::solve_S()
{
    structuralMatrix();

    for (size_t j = 1; j < ny-1; j++)
    {
        // BC west
        zBoundary(z.westBC, z.west(j),    0, j, z.r1West, z.r2West, h_1s1_west(j), h_1s2_west(j));
        // BC east
        zBoundary(z.eastBC, z.east(j), nx-1, j, z.r1East, z.r2East, h_1s1_east(j), h_1s2_east(j));
    }
    for (size_t i = 1; i < nx-1; i++)
    {
        // BC south
        zBoundary(z.southBC, z.south(i), i,    0, z.r1South, z.r2South, h_2s1_south(i), h_2s2_south(i));
        // BC north
        zBoundary(z.northBC, z.north(i), i, ny-1, z.r1North, z.r2North, h_2s1_north(i), h_2s2_north(i));
    }
    // BC south-west corner (i = 0, j = 0)
    if (chi[0]->curveType == CurveType::Interface && chi[3]->curveType == CurveType::Boundary)
        zBoundary(z.westBC, z.west(0), 0, 0, z.r1West, z.r2West, h_1s1_west(0), h_1s2_west(0));
    else if (chi[0]->curveType == CurveType::Boundary  && chi[3]->curveType == CurveType::Interface)
        zBoundary(z.southBC, z.south(0), 0, 0, z.r1South, z.r2South, h_2s1_south(0), h_2s2_south(0));
    else
        zBoundary(z.southBC, z.south(0), 0, 0, z.r1South, z.r2South, h_2s1_south(0), h_2s2_south(0));

    // BC north-west corner (i = 0, j = ny-1)
    size_t j = ny-1;
    size_t k = j*nx;
    if (chi[2]->curveType == CurveType::Interface && chi[3]->curveType == CurveType::Boundary)
        zBoundary(z.westBC, z.west(j), 0, j, z.r1West, z.r2West, h_1s1_west(j), h_1s2_west(j));
    else if (chi[2]->curveType == CurveType::Boundary  && chi[3]->curveType == CurveType::Interface)
        zBoundary(z.northBC, z.north(0), 0, j, z.r1North, z.r2North, h_2s1_north(0), h_2s2_north(0));
    else
        zBoundary(z.westBC,  z.west(j),  0, j, z.r1West,  z.r2West,  h_1s1_west(j),  h_1s2_west(j));

    // BC south-east corner (i = nx-1, j = 0)
    size_t i = nx-1;
    if (chi[0]->curveType == CurveType::Interface && chi[1]->curveType == CurveType::Boundary)
        zBoundary(z.eastBC, z.east(0), i, 0, z.r1East, z.r2East, h_1s1_east(0), h_1s2_east(0));
    else if (chi[0]->curveType == CurveType::Boundary  && chi[1]->curveType == CurveType::Interface)
        zBoundary(z.southBC, z.south(i), i, 0, z.r1South, z.r2South, h_2s1_south(i), h_2s2_south(i));
    else
        zBoundary(z.eastBC,  z.east(0),  i, 0, z.r1East,  z.r2East,  h_1s1_east(0),  h_1s2_east(0));

    // BC north-east corner (i = nx-1, j = ny-1)
    i = nx-1;
    j = ny-1;
    k = i + j*nx;
    if (chi[2]->curveType == CurveType::Interface && chi[1]->curveType == CurveType::Boundary)
        zBoundary(z.eastBC, z.east(j), i, j, z.r1East, z.r2East, h_1s1_east(j), h_1s2_east(j));
    else if (chi[2]->curveType == CurveType::Boundary  && chi[1]->curveType == CurveType::Interface)
        zBoundary(z.northBC, z.north(i), i, j, z.r1North, z.r2North, h_2s1_north(i), h_2s2_north(i));
    else
        zBoundary(z.northBC, z.north(i), i, j, z.r1North, z.r2North, h_2s1_north(i), h_2s2_north(i));

    arma::lu(L, U, P, S);
}

void Membrane::solve_b()
{
    #pragma omp parallel for
    for (size_t j = 1; j < ny-1; j++)
        for (size_t i = 1; i < nx-1; i++)
            b(j*nx+i) =-p(i, j);

    z.arma::mat::operator=(reshape(solve(trimatu(U), solve(trimatl(L), P*b)), nx, ny));
}