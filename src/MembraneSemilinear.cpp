#include "Membrane.hpp"

void Membrane::semilinear()
{
    analysis = Analysis::semilinear;
    structuralMatrix();
    arma::mat e11  = ec.slice(0);
    arma::mat e12  = ec.slice(1);
    arma::mat e22  = ec.slice(2);
    arma::vec E11  = vectorise(e11);
    arma::vec E12  = vectorise(e12);
    arma::vec E22  = vectorise(e22);

    for (size_t substep = 1; substep <= substeps; substep++)
    {
        printf("Substep %lu/%lu\n", substep, substeps);
        arma::mat pressure = p*substep/substeps;
        for (size_t q = 0; q < iter; q++)
        {
            printf("Iteration %lu/%lu\n", q+1, iter);
            arma::vec Z  = vectorise(arma::mat(z));
            arma::vec P  = vectorise(pressure);
            double pMAX = max(arma::abs(P));
            double pRMS = sqrt(sum(P%P));
            arma::vec z_1 = DD1*Z;
            arma::vec z_2 = DD2*Z;
            arma::vec sae = sqrt(1 + E11%z_1%z_1 + 2*E12%z_1%z_2 + E22%z_2%z_2);
            b = S*Z + P%sae;
            arma::mat dsaedz = (QLM((E11%z_1 + E12%z_2), DD1) + QLM((E12%z_1 + E22%z_2), DD2))&(sae);

            arma::mat A = S + QLM(P, dsaedz);
            for (size_t j = 1; j < ny-1; j++)
            {
                zBoundarySemilinear(z.westBC, z.west(j),    0, j, A, b, z.r1West, z.r2West, z_1(j*nx),       z_2(j*nx),       h_1s1_west(j), h_1s2_west(j));
                zBoundarySemilinear(z.eastBC, z.east(j), nx-1, j, A, b, z.r1East, z.r2East, z_1((j+1)*nx-1), z_2((j+1)*nx-1), h_1s1_east(j), h_1s2_east(j));
            }
            for (size_t i = 1; i < nx-1; i++)
            {
                zBoundarySemilinear(z.southBC, z.south(i), i,    0, A, b, z.r1South, z.r2South, z_1(i),           z_2(i),           h_2s1_south(i), h_2s2_south(i));
                zBoundarySemilinear(z.northBC, z.north(i), i, ny-1, A, b, z.r1North, z.r2North, z_1(i+(ny-1)*nx), z_2(i+(ny-1)*nx), h_2s1_north(i), h_2s2_north(i));
            }
            // BC south-west corner (i = 0, j = 0)
            zBoundarySemilinear(z.southBC, z.south(0), 0, 0, A, b, z.r1South, z.r2South, z_1(0), z_2(0), h_2s1_south(0), h_2s2_south(0));

            // BC north-west corner (i = 0, j = ny-1)
            size_t j = ny-1;
            size_t k = j*nx;
            zBoundarySemilinear(z.westBC,  z.west(j),  0, j, A, b, z.r1West,  z.r2West,  z_1(k), z_2(k), h_1s1_west(j),  h_1s2_west(j));

            // BC south-east corner (i = nx-1, j = 0)
            size_t i = nx-1;
            zBoundarySemilinear(z.eastBC,  z.east(0),  i, 0, A, b, z.r1East,  z.r2East,  z_1(i), z_2(i), h_1s1_east(0),  h_1s2_east(0));

            // BC north-east corner (i = nx-1, j = ny-1)
            i = nx-1;
            j = ny-1;
            k = i + j*nx;
            zBoundarySemilinear(z.northBC, z.north(i), i, j, A, b, z.r1North, z.r2North, z_1(k), z_2(k), h_2s1_north(i), h_2s2_north(i));

            double resMAX = max(arma::abs(b))/pMAX;
            double resRMS = sqrt(sum(b%b)/nxy)/pRMS;
            std::println("MAX Residual: {0:4.2e}",   resMAX);
            std::println("RMS Residual: {0:4.2e}\n", resRMS);
            if (resRMS < residualTarget)
                break;

            arma::vec dz = solve(A, b);
            double omega = armijoSemilinear(dz, P);
            std::println("omega = {0:4.2f}", omega);
            z.arma::mat::operator-=(omega*reshape(dz, nx, ny));
        }
    }
}

double Membrane::armijoSemilinear(arma::vec dz, arma::vec &P)
{
    double omega = 1;
    double T = residualLevelFunctionSemilinear(vectorise(arma::mat(z)), P);
    double T_dx_min = 1e10;
    for (size_t m = 0; m < 10; m++)
    {
        double alpha = pow(0.5, m);
        double T_dx = residualLevelFunctionSemilinear(vectorise(arma::mat(z)) - alpha*dz, P);
        if (T_dx <= (1 - alpha/2)*T)
            if (T_dx < T_dx_min)
            {
                omega = alpha;
                T_dx_min = T_dx;
            }
    }
    return omega;
}

double Membrane::residualLevelFunctionSemilinear(arma::vec Z, arma::vec &P)
{
    arma::vec Z_1 = DD1*Z;
    arma::vec Z_2 = DD2*Z;
    arma::vec f = S*Z + P%(1 + vectorise(e11())%Z_1%Z_1 + 2*vectorise(e12())%Z_1%Z_2 + vectorise(e22())%Z_2%Z_2);

    return dot(f, f);
}