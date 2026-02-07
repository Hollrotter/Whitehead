#include "Membrane.hpp"

void Membrane::planeStrain()
{
    size_t nxy2 = 2*nxy;

    arma::mat Phi111 = 3*H1111%g111() + 6*H1112%g112() + (H1122 + 2*H1221)%g122();
    arma::mat Phi112 = 2*H1112%(g111() + g212()) + 2*H1222%g122() + H1111%g211() + H1221%(3*g112() + g222()) + H1122%g112();
    arma::mat Phi121 = 2*H1111%g211() + (H1122 + H1221)%(g112() + g222()) + H1112%(g111() + 4*g212()) + H1222%g122();
    arma::mat Phi122 = 3*H1112%g211() + 2*(2*H1221 + H1122)%g212() + 3*H1222%g222();
    arma::mat Psi11  = H1111%Psi1111 + H1221%Psi1221 + (H1221 + H1122)%Psi2121 + H1222%Psi2221 + H1112%(Psi2112 + 2*Psi1121);
    arma::mat Psi12  = H1111%Psi1112 + H1221%Psi1222 + H1112%(Psi2111 + 2*Psi1122) + (H1221 + H1122)%Psi2122 + H1222%Psi2222;

    arma::mat Phi211 = 3*H1112%g111() + 3*H1222%g122() + 2*(2*H1221 + H1122)%g112();
    arma::mat Phi212 = (H1122 + H1221)%(g111() + g212()) + H1112%g211() + H1222%(4*g112() + g222()) + 2*H2222%g122();
    arma::mat Phi221 = 2*H1112%g211() + H1122%g212() + 2*H1222%(g112() + g222()) + H1221%(g111() + 3*g212()) + H2222%g122();
    arma::mat Phi222 = 6*H1222%g212() + (H1122 + 2*H1221)%g211() + 3*H2222%g222();
    arma::mat Psi21  = H1112%Psi1111 + H1221%Psi2112 + H1222%(Psi1221 + 2*Psi2121) + (H1122 + H1221)%Psi1121 + H2222%Psi2221;
    arma::mat Psi22  = H1112%Psi1112 + H1221%Psi2111 + (H1221 + H1122)%Psi1122 + H2222%Psi2222 + H1222%(Psi1222 + 2*Psi2122);

    arma::mat A = join_vert(join_horiz(d2dx2(H1111) + d2dxdy(2*H1112)       + d2dy2(H1221) - ddx(Phi111) - ddy(Phi112) + constant(Psi11),
                                       d2dx2(H1112) + d2dxdy(H1122 + H1221) + d2dy2(H1222) - ddx(Phi121) - ddy(Phi122) + constant(Psi12)),
                            join_horiz(d2dx2(H1112) + d2dxdy(H1122 + H1221) + d2dy2(H1222) - ddx(Phi211) - ddy(Phi212) + constant(Psi21),
                                       d2dx2(H1221) + d2dxdy(2*H1222)       + d2dy2(H2222) - ddx(Phi221) - ddy(Phi222) + constant(Psi22)));
                                       
    arma::vec bv = join_vert(-p1,-p2)/D;
    
    for (size_t j = 1; j < ny-1; j++)
    {
        // BC west

        v1BoundaryWestEast(v1.westBC, v1.west(j), 0, j, A, bv, v1.r1West, v1.r2West, h_1s1_west(j), h_1s2_west(j));
        v2BoundaryWestEast(v2.westBC, v2.west(j), 0, j, A, bv, v2.r1West, v2.r2West, h_1s1_west(j), h_1s2_west(j), h_2s2_west(j));
        n11BoundaryLinear(n11.westBC, n11.west(j), 0, j, A, bv, h_11s_west(j));
        n12BoundaryLinear(n12.westBC, n12.west(j), 0, j, A, bv, h_11s_west(j), h_12s_west(j), h_22s_west(j));

        // BC east
        size_t i = nx-1;

        v1BoundaryWestEast(v1.eastBC, v1.east(j), i, j, A, bv, v1.r1East, v1.r2East, h_1s1_east(j), h_1s2_east(j));
        v2BoundaryWestEast(v2.eastBC, v2.east(j), i, j, A, bv, v2.r1East, v2.r2East, h_1s1_east(j), h_1s2_east(j), h_2s2_east(j));
        n11BoundaryLinear(n11.eastBC, n11.east(j), i, j, A, bv, h_11s_east(j));
        n12BoundaryLinear(n12.eastBC, n12.east(j), i, j, A, bv, h_11s_east(j), h_12s_east(j), h_22s_east(j));
    }

    for (size_t i = 1; i < nx-1; i++)
    {
        // BC south

        v1BoundarySouthNorth(v1.southBC, v1.south(i), i, 0, A, bv, v1.r1South, v1.r2South, h_1s1_south(i), h_2s1_south(i), h_2s2_south(i));
        v2BoundarySouthNorth(v2.southBC, v2.south(i), i, 0, A, bv, v2.r1South, v2.r2South, h_2s1_south(i), h_2s2_south(i));
        n21BoundaryLinear(n12.southBC, n12.south(i), i, 0, A, bv, h_11s_south(i), h_21s_south(i), h_22s_south(i));
        n22BoundaryLinear(n22.southBC, n22.south(i), i, 0, A, bv, h_22s_south(i));
        
        // BC north
        size_t j = ny-1;

        v1BoundarySouthNorth(v1.northBC, v1.north(i), i, j, A, bv, v1.r1North, v1.r2North, h_1s1_north(i), h_2s1_north(i), h_2s2_north(i));
        v2BoundarySouthNorth(v2.northBC, v2.north(i), i, j, A, bv, v2.r1North, v2.r2North, h_2s1_north(i), h_2s2_north(i));
        n21BoundaryLinear(n12.northBC, n12.north(i), i, j, A, bv, h_11s_north(i), h_21s_north(i), h_22s_north(i));
        n22BoundaryLinear(n22.northBC, n22.north(i), i, j, A, bv, h_22s_north(i));
    }

    // BC south-west corner (i = 0, j = 0)
    v1BoundarySouthNorth(v1.southBC, v1.south(0), 0, 0, A, bv, v1.r1South, v1.r2South, h_1s1_south(0), h_2s1_south(0), h_2s2_south(0));
    v2BoundarySouthNorth(v2.southBC, v2.south(0), 0, 0, A, bv, v2.r1South, v2.r2South, h_2s1_south(0), h_2s2_south(0));
    n21BoundaryLinear(n12.southBC, n12.south(0), 0, 0, A, bv, h_11s_south(0), h_21s_south(0), h_22s_south(0));
    n22BoundaryLinear(n22.southBC, n22.south(0), 0, 0, A, bv, h_22s_south(0));

    // BC north-west corner (i = 0, j = ny-1)
    size_t j = ny-1;
    size_t k = j*nx;
    v1BoundaryWestEast(v1.westBC, v1.west(j), 0, j, A, bv, v1.r1West, v1.r2West, h_1s1_west(j), h_1s2_west(j));
    v2BoundaryWestEast(v2.westBC, v2.west(j), 0, j, A, bv, v2.r1West, v2.r2West, h_1s1_west(j), h_1s2_west(j), h_2s2_west(j));
    n11BoundaryLinear(n11.westBC, n11.west(j), 0, j, A, bv, h_11s_west(j));
    n12BoundaryLinear(n12.westBC, n12.west(j), 0, j, A, bv, h_11s_west(j), h_12s_west(j), h_22s_west(j));

    // BC south-east corner (i = nx-1, j = 0)
    size_t i = nx-1;
    v1BoundaryWestEast(v1.eastBC, v1.east(0), i, 0, A, bv, v1.r1East, v1.r2East, h_1s1_east(0), h_1s2_east(0));
    v2BoundaryWestEast(v2.eastBC, v2.east(0), i, 0, A, bv, v2.r1East, v2.r2East, h_1s1_east(0), h_1s2_east(0), h_2s2_east(0));
    n11BoundaryLinear(n11.eastBC, n11.east(0), i, 0, A, bv, h_11s_east(0));
    n12BoundaryLinear(n12.eastBC, n12.east(0), i, 0, A, bv, h_11s_east(0), h_12s_east(0), h_22s_east(0));

    // BC north-east corner (i = nx-1, j = ny-1)
    i = nx-1;
    j = ny-1;
    k = i + j*nx;
    v1BoundarySouthNorth(v1.northBC, v1.north(i), i, j, A, bv, v1.r1North, v1.r2North, h_1s1_north(i), h_2s1_north(i), h_2s2_north(i));
    v2BoundarySouthNorth(v2.northBC, v2.north(i), i, j, A, bv, v2.r1North, v2.r2North, h_2s1_north(i), h_2s2_north(i));
    n21BoundaryLinear(n12.northBC, n12.north(i), i, j, A, bv, h_11s_north(i), h_21s_north(i), h_22s_north(i));
    n22BoundaryLinear(n22.northBC, n22.north(i), i, j, A, bv, h_22s_north(i));

    bv = solve(A, bv);
    v1 = arma::reshape(bv(arma::span(  0,   nxy-1)), nx, ny);
    v2 = arma::reshape(bv(arma::span(nxy, 2*nxy-1)), nx, ny);

    v1__1 = D1*v1     - g111()%v1 - g211()%v2;
    v1__2 = v1*D2.t() - g112()%v1 - g212()%v2;
    v2__1 = D1*v2     - g112()%v1 - g212()%v2;
    v2__2 = v2*D2.t() - g122()%v1 - g222()%v2;

    gamma_11 = v1__1;
    gamma_12 = (v1__2 + v2__1)/2;
    gamma_22 = v2__2;

    n11 = D*(H1111%gamma_11 + 2*H1112%gamma_12 + H1122%gamma_22);
    n12 = D*(H1112%gamma_11 + 2*H1221%gamma_12 + H1222%gamma_22);
    n22 = D*(H1122%gamma_11 + 2*H1222%gamma_12 + H2222%gamma_22);
}