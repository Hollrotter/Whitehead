#include "Membrane.hpp"

void Membrane::nonlinear()
{
    analysis = Analysis::nonlinear;

    arma::mat D2_t  = D2.t();
    arma::mat D22_t = D22.t();

    arma::mat G112DD1 = QLM(2*g112(), DD1);
    arma::mat G212DD2 = QLM(2*g212(), DD2);
    arma::mat G211DD2 = QLM(  g211(), DD2);
    arma::mat G122DD1 = QLM(  g122(), DD1);

    arma::mat dz__11dz = D2D11 - QLM(g111(), DD1) -     G211DD2;
    arma::mat dz__12dz = D2D12 - QLM(g112(), DD1) - QLM(g212(), DD2);
    arma::mat dz__22dz = D2D22 -     G122DD1  - QLM(g222(), DD2);

    arma::mat dv1__1dv1  = DD1 - diagmat(vectorise(g111()));
    arma::mat dv1__2dv1  = DD2 - diagmat(vectorise(g112()));
    arma::mat dv1__1dv2  =     - diagmat(vectorise(g211()));
    arma::mat dv1__2dv2  =     - diagmat(vectorise(g212()));
    arma::mat dv2__1dv1  =     - diagmat(vectorise(g112()));
    arma::mat dv2__2dv1  =     - diagmat(vectorise(g122()));
    arma::mat dv2__1dv2  = DD1 - diagmat(vectorise(g212()));
    arma::mat dv2__2dv2  = DD2 - diagmat(vectorise(g222()));

    arma::mat dv1__11dv1 = D2D11 - QLM(3*g111(), DD1) - G211DD2 + diagmat(vectorise(Psi1111));
    arma::mat dv2__22dv2 = D2D22 - QLM(3*g222(), DD2) - G122DD1 + diagmat(vectorise(Psi2222));

    arma::mat dv1__11dv2 = diagmat(vectorise(Psi1112)) - QLM(2*g211(), DD1);
    arma::mat dv2__22dv1 = diagmat(vectorise(Psi2221)) - QLM(2*g122(), DD2);

    arma::mat dv1__12dv1 = D2D12 - G112DD1 - QLM((g111() + g212()), DD2) + diagmat(vectorise(Psi1121));
    arma::mat dv2__12dv2 = D2D12 - G212DD2 - QLM((g112() + g222()), DD1) + diagmat(vectorise(Psi2122));

    arma::mat dv1__12dv2 = diagmat(vectorise(Psi1122)) - G211DD2 - QLM(g212(), DD1);
    arma::mat dv2__12dv1 = diagmat(vectorise(Psi2121)) - G122DD1 - QLM(g112(), DD2);

    arma::mat dv1__22dv1 = D2D22 - G122DD1 - QLM((2*g112() - g222()), DD2) + diagmat(vectorise(Psi1221));
    arma::mat dv2__11dv2 = D2D11 - G211DD2 - QLM((2*g212() + g111()), DD1) + diagmat(vectorise(Psi2111));

    arma::mat dv1__22dv2 = diagmat(vectorise(Psi1222)) - G212DD2;
    arma::mat dv2__11dv1 = diagmat(vectorise(Psi2112)) - G112DD1;

    size_t nxy2 = 2*nxy;
    size_t nxy3 = 3*nxy;

    arma::mat v1_1 = D1*v1;
    arma::mat v1_2 = v1*D2_t;
    arma::mat v2_1 = D1*v2;
    arma::mat v2_2 = v2*D2_t;
    v1__1 = v1_1 - g111()%v1 - g211()%v2;
    v1__2 = v1_2 - g112()%v1 - g212()%v2;
    v2__1 = v2_1 - g112()%v1 - g212()%v2;
    v2__2 = v2_2 - g122()%v1 - g222()%v2;
    arma::mat z_1  = D1*z;
    arma::mat z_2  = z*D2_t;

    gamma_11 =  v1__1         + (z_1%z_1 + e11()%v1__1%v1__1 + 2*e12()%v1__1%v2__1                 + e22()%v2__1%v2__1)/2;
    gamma_12 = (v1__2 + v2__1 +  z_1%z_2 + e11()%v1__1%v1__2 +   e12()%(v1__1%v2__2 + v2__1%v1__2) + e22()%v2__1%v2__2)/2;
    gamma_22 =  v2__2         + (z_2%z_2 + e11()%v1__2%v1__2 + 2*e12()%v1__2%v2__2                 + e22()%v2__2%v2__2)/2;

    n11 = D*(H1111%gamma_11 + 2*H1112%gamma_12 + H1122%gamma_22);
    n12 = D*(H1112%gamma_11 + 2*H1221%gamma_12 + H1222%gamma_22);
    n22 = D*(H1122%gamma_11 + 2*H1222%gamma_12 + H2222%gamma_22);

    for (size_t substep = 1; substep <= substeps; substep++)
    {
        printf("Substep %lu/%lu\n", substep, substeps);
        arma::mat pressure =  p*substep/substeps;
        arma::mat surface1 = reshape(p1, nx, ny)*substep/substeps;
        arma::mat surface2 = reshape(p2, nx, ny)*substep/substeps;
        for (size_t q = 0; q < iter; q++)
        {
            printf("Iteration %lu/%lu\n", q+1, iter);
            double pMAX = max(max(arma::abs(pressure)));
            double pRMS = sqrt(sum(sum(pressure%pressure)));

            arma::mat z_11 = D11*z;
            arma::mat z_12 = D1*z*D2_t;
            arma::mat z_22 = z*D22_t;

            arma::mat v1_11 = D11*v1;
            arma::mat v1_12 = D1*v1*D2_t;
            arma::mat v1_22 = v1*D22_t;
            arma::mat v2_11 = D11*v2;
            arma::mat v2_12 = D1*v2*D2_t;
            arma::mat v2_22 = v2*D22_t;

            arma::mat v1__11 = v1_11 - 3*g111()%v1_1 - g211()%(v1_2 + 2*v2_1) + Psi1111%v1 + Psi1112%v2;
            arma::mat v2__22 = v2_22 - 3*g222()%v2_2 - g122()%(v2_1 + 2*v1_2) + Psi2221%v1 + Psi2222%v2;
            arma::mat v1__12 = v1_12 - 2*g112()%v1_1 - (g111() + g212())%v1_2 - g212()%v2_1 - g211()%v2_2 + Psi1121%v1 + Psi1122%v2;
            arma::mat v2__12 = v2_12 - 2*g212()%v2_2 - (g112() + g222())%v2_1 - g122()%v1_1 - g112()%v1_2 + Psi2121%v1 + Psi2122%v2;
            arma::mat v1__22 = v1_22 -   g122()%v1_1 - (2*g112() + g222())%v1_2 - 2*g212()%v2_2 + Psi1221%v1 + Psi1222%v2;
            arma::mat v2__11 = v2_11 -   g211()%v2_2 - (2*g212() + g111())%v2_1 - 2*g112()%v1_1 + Psi2112%v1 + Psi2111%v2;

            arma::mat z__11 = z_11 - g111()%z_1 - g211()%z_2;
            arma::mat z__12 = z_12 - g112()%z_1 - g212()%z_2;
            arma::mat z__22 = z_22 - g122()%z_1 - g222()%z_2;

            arma::mat z_1__12 = z_1%z__12;
            arma::mat z_2__12 = z_2%z__12;

            arma::mat z_112 = z_1__12 + z_2%z__11;
            arma::mat z_122 = z_1%z__22 + z_2__12;

            arma::mat G_11__1 =  v1__11        + z_1%z__11 + e11()%v1__1%v1__11                  + e12()%(v1__1%v2__11 + v2__1%v1__11)                               + e22()%v2__1%v2__11;
            arma::mat G_11__2 =  v1__12        + z_1__12   + e11()%v1__1%v1__12                  + e12()%(v1__1%v2__12 + v2__1%v1__12)                               + e22()%v2__1%v2__12;
            arma::mat G_12__1 = (v1__12 + v2__11 + z_112   + e11()%(v1__1%v1__12 + v1__2%v1__11) + e12()%(v2__1%v1__12 + v1__2%v2__11 + v1__1%v2__12 + v2__2%v1__11) + e22()%(v2__1%v2__12 + v2__2%v2__11))/2;
            arma::mat G_12__2 = (v1__22 + v2__12 + z_122   + e11()%(v1__1%v1__22 + v1__2%v1__12) + e12()%(v2__1%v1__22 + v1__2%v2__12 + v1__1%v2__22 + v2__2%v1__12) + e22()%(v2__1%v2__22 + v2__2%v2__12))/2;
            arma::mat G_22__1 =  v2__12        + z_2__12   + e11()%v1__2%v1__12                  + e12()%(v1__2%v2__12 + v2__2%v1__12)                               + e22()%v2__2%v2__12;
            arma::mat G_22__2 =  v2__22        + z_2%z__22 + e11()%v1__2%v1__22                  + e12()%(v1__2%v2__22 + v2__2%v1__22)                               + e22()%v2__2%v2__22;

            arma::mat n11__1 = D*(H1111%G_11__1 + 2*H1112%G_12__1 + H1122%G_22__1);
            arma::mat n12__1 = D*(H1112%G_11__1 + 2*H1221%G_12__1 + H1222%G_22__1);
            arma::mat n12__2 = D*(H1112%G_11__2 + 2*H1221%G_12__2 + H1222%G_22__2);
            arma::mat n22__2 = D*(H1122%G_11__2 + 2*H1222%G_12__2 + H2222%G_22__2);

            arma::mat L111 = e11()%v1__11 + e12()%v2__11;
            arma::mat L112 = e11()%v1__12 + e12()%v2__12;
            arma::mat L122 = e11()%v1__22 + e12()%v2__22;
            arma::mat L211 = e12()%v1__11 + e22()%v2__11;
            arma::mat L212 = e12()%v1__12 + e22()%v2__12;
            arma::mat L222 = e12()%v1__22 + e22()%v2__22;

            arma::mat w   =-(e11()%z_1%z_1 + 2*e12()%z_1%z_2 + e22()%z_2%z_2)/2;
            arma::mat w_1 = z_1%(e11()%v1__1 + e12()%v2__1 - 1) + z_2%(e12()%v1__1 + e22()%v2__1);
            arma::mat w_2 = z_1%(e11()%v1__2 + e12()%v2__2)     + z_2%(e12()%v1__2 + e22()%v2__2 - 1);

            arma::mat w1  = e11()%w_1 + e12()%w_2;
            arma::mat w2  = e12()%w_1 + e22()%w_2;
            
            arma::mat s1 = surface1%(e11()%v1__1 + e12()%v2__1 + 1) + surface2%(e11()%v1__2 + e12()%v2__2)     + pressure%w1;
            arma::mat s2 = surface1%(e12()%v1__1 + e22()%v2__1)     + surface2%(e12()%v1__2 + e22()%v2__2 + 1) + pressure%w2;
            arma::mat s3 = surface1%z_1                             + surface2%z_2                             + pressure%(1 + w);

            arma::mat sae = 1 + e11()%gamma_11 + 2*e12()%gamma_12 + e22()%gamma_22;
            b = join_vert(vectorise(n11__1 + n12__2 + L111%n11 + 2*L112%n12 + L122%n22 + sae%s1)/D,
                          vectorise(n12__1 + n22__2 + L211%n11 + 2*L212%n12 + L222%n22 + sae%s2)/D,
                          vectorise(z__11%n11 + 2*z__12%n12 + z__22%n22 + sae%s3)/D);
            
            arma::mat dz_1__12 = QLM(z_1, dz__12dz);
            arma::mat dz_2__12 = QLM(z_2, dz__12dz);

            arma::mat dG_11__1dz  = QLM(z__11, DD1) + QLM(z_1, dz__11dz);
            arma::mat dG_11__2dz  = QLM(z__12, DD1) + dz_1__12;
            arma::mat dG_12__1dz  = (QLM(z__12, DD1) + dz_1__12 + QLM(z__11, DD2) + QLM(z_2, dz__11dz))/2;
            arma::mat dG_12__2dz  = (QLM(z__22, DD1) + QLM(z_1, dz__22dz) + QLM(z__12, DD2) + dz_2__12)/2;
            arma::mat dG_22__1dz  = QLM(z__12, DD2) + dz_2__12;
            arma::mat dG_22__2dz  = QLM(z__22, DD2) + QLM(z_2, dz__22dz);

            arma::mat dG_11dv1 = dv1__1dv1 + QLM(e11()%v1__1, dv1__1dv1) + QLM(e12()%v1__1, dv2__1dv1) + QLM(e12()%v2__1, dv1__1dv1) + QLM(e22()%v2__1, dv2__1dv1);
            arma::mat dG_11dv2 = dv1__1dv2 + QLM(e11()%v1__1, dv1__1dv2) + QLM(e12()%v1__1, dv2__1dv2) + QLM(e12()%v2__1, dv1__1dv2) + QLM(e22()%v2__1, dv2__1dv2);
            arma::mat dG_12dv1 = (dv1__2dv1 + dv2__1dv1 + QLM(e11()%v1__1, dv1__2dv1) + QLM(e11()%v1__2, dv1__1dv1) + QLM(e12()%v1__1, dv2__2dv1) + QLM(e12()%v2__2, dv1__1dv1) + QLM(e12()%v1__2, dv2__1dv1) + QLM(e12()%v2__1, dv1__2dv1) + QLM(e22()%v2__1, dv2__2dv1) + QLM(e22()%v2__2, dv2__1dv1))/2;
            arma::mat dG_12dv2 = (dv1__2dv2 + dv2__1dv2 + QLM(e11()%v1__1, dv1__2dv2) + QLM(e11()%v1__2, dv1__1dv2) + QLM(e12()%v1__1, dv2__2dv2) + QLM(e12()%v2__2, dv1__1dv2) + QLM(e12()%v1__2, dv2__1dv2) + QLM(e12()%v2__1, dv1__2dv2) + QLM(e22()%v2__1, dv2__2dv2) + QLM(e22()%v2__2, dv2__1dv2))/2; 
            arma::mat dG_22dv1 = dv2__2dv1 + QLM(e11()%v1__2, dv1__2dv1) + QLM(e12()%v1__2, dv2__2dv1) + QLM(e12()%v2__2, dv1__2dv1) + QLM(e22()%v2__2, dv2__2dv1);
            arma::mat dG_22dv2 = dv2__2dv2 + QLM(e11()%v1__2, dv1__2dv2) + QLM(e12()%v1__2, dv2__2dv2) + QLM(e12()%v2__2, dv1__2dv2) + QLM(e22()%v2__2, dv2__2dv2);

            arma::mat dG_11__1dv1 = dv1__11dv1 + QLM(e11()%v1__1, dv1__11dv1) + QLM(e11()%v1__11, dv1__1dv1) + QLM(e12()%v2__1, dv1__11dv1) + QLM(e12()%v1__11, dv2__1dv1) + QLM(e12()%v1__1, dv2__11dv1) + QLM(e12()%v2__11, dv1__1dv1) + QLM(e22()%v2__1, dv2__11dv1) + QLM(e22()%v2__11, dv2__1dv1);
            arma::mat dG_11__1dv2 = dv1__11dv2 + QLM(e11()%v1__1, dv1__11dv2) + QLM(e11()%v1__11, dv1__1dv2) + QLM(e12()%v2__1, dv1__11dv2) + QLM(e12()%v1__11, dv2__1dv2) + QLM(e12()%v1__1, dv2__11dv2) + QLM(e12()%v2__11, dv1__1dv2) + QLM(e22()%v2__1, dv2__11dv2) + QLM(e22()%v2__11, dv2__1dv2);
            arma::mat dG_11__2dv1 = dv1__12dv1 + QLM(e11()%v1__1, dv1__12dv1) + QLM(e11()%v1__12, dv1__1dv1) + QLM(e12()%v2__1, dv1__12dv1) + QLM(e12()%v1__12, dv2__1dv1) + QLM(e12()%v1__1, dv2__12dv1) + QLM(e12()%v2__12, dv1__1dv1) + QLM(e22()%v2__1, dv2__12dv1) + QLM(e22()%v2__12, dv2__1dv1);
            arma::mat dG_11__2dv2 = dv1__12dv2 + QLM(e11()%v1__1, dv1__12dv2) + QLM(e11()%v1__12, dv1__1dv2) + QLM(e12()%v2__1, dv1__12dv2) + QLM(e12()%v1__12, dv2__1dv2) + QLM(e12()%v1__1, dv2__12dv2) + QLM(e12()%v2__12, dv1__1dv2) + QLM(e22()%v2__1, dv2__12dv2) + QLM(e22()%v2__12, dv2__1dv2);

            arma::mat dG_12__1dv1 = (dv1__12dv1 + dv2__11dv1 + QLM(e11()%v1__1, dv1__12dv1) + QLM(e11()%v1__12, dv1__1dv1) + QLM(e11()%v1__2, dv1__11dv1) + QLM(e11()%v1__11, dv1__2dv1) + QLM(e12()%v2__1, dv1__12dv1) + QLM(e12()%v1__12, dv2__1dv1) + QLM(e12()%v1__2, dv2__11dv1) + QLM(e12()%v2__11, dv1__2dv1) + QLM(e12()%v1__1, dv2__12dv1) + QLM(e12()%v2__12, dv1__1dv1) + QLM(e12()%v2__2, dv1__11dv1) + QLM(e12()%v1__11, dv2__2dv1) + QLM(e22()%v2__1, dv2__12dv1) + QLM(e22()%v2__12, dv2__1dv1) + QLM(e22()%v2__2, dv2__11dv1) + QLM(e22()%v2__11, dv2__2dv1))/2;
            arma::mat dG_12__1dv2 = (dv1__12dv2 + dv2__11dv2 + QLM(e11()%v1__1, dv1__12dv2) + QLM(e11()%v1__12, dv1__1dv2) + QLM(e11()%v1__2, dv1__11dv2) + QLM(e11()%v1__11, dv1__2dv2) + QLM(e12()%v2__1, dv1__12dv2) + QLM(e12()%v1__12, dv2__1dv2) + QLM(e12()%v1__2, dv2__11dv2) + QLM(e12()%v2__11, dv1__2dv2) + QLM(e12()%v1__1, dv2__12dv2) + QLM(e12()%v2__12, dv1__1dv2) + QLM(e12()%v2__2, dv1__11dv2) + QLM(e12()%v1__11, dv2__2dv2) + QLM(e22()%v2__1, dv2__12dv2) + QLM(e22()%v2__12, dv2__1dv2) + QLM(e22()%v2__2, dv2__11dv2) + QLM(e22()%v2__11, dv2__2dv2))/2;
            arma::mat dG_12__2dv1 = (dv1__22dv1 + dv2__12dv1 + QLM(e11()%v1__1, dv1__22dv1) + QLM(e11()%v1__22, dv1__1dv1) + QLM(e11()%v1__2, dv1__12dv1) + QLM(e11()%v1__12, dv1__2dv1) + QLM(e12()%v2__1, dv1__22dv1) + QLM(e12()%v1__22, dv2__1dv1) + QLM(e12()%v1__2, dv2__12dv1) + QLM(e12()%v2__12, dv1__2dv1) + QLM(e12()%v1__1, dv2__22dv1) + QLM(e12()%v2__22, dv1__1dv1) + QLM(e12()%v2__2, dv1__12dv1) + QLM(e12()%v1__12, dv2__2dv1) + QLM(e22()%v2__1, dv2__22dv1) + QLM(e22()%v2__22, dv2__1dv1) + QLM(e22()%v2__2, dv2__12dv1) + QLM(e22()%v2__12, dv2__2dv1))/2;
            arma::mat dG_12__2dv2 = (dv1__22dv2 + dv2__12dv2 + QLM(e11()%v1__1, dv1__22dv2) + QLM(e11()%v1__22, dv1__1dv2) + QLM(e11()%v1__2, dv1__12dv2) + QLM(e11()%v1__12, dv1__2dv2) + QLM(e12()%v2__1, dv1__22dv2) + QLM(e12()%v1__22, dv2__1dv2) + QLM(e12()%v1__2, dv2__12dv2) + QLM(e12()%v2__12, dv1__2dv2) + QLM(e12()%v1__1, dv2__22dv2) + QLM(e12()%v2__22, dv1__1dv2) + QLM(e12()%v2__2, dv1__12dv2) + QLM(e12()%v1__12, dv2__2dv2) + QLM(e22()%v2__1, dv2__22dv2) + QLM(e22()%v2__22, dv2__1dv2) + QLM(e22()%v2__2, dv2__12dv2) + QLM(e22()%v2__12, dv2__2dv2))/2;

            arma::mat dG_22__1dv1 = dv2__12dv1 + QLM(e11()%v1__2, dv1__12dv1) + QLM(e11()%v1__12, dv1__2dv1) + QLM(e12()%v2__2, dv1__12dv1) + QLM(e12()%v1__12, dv2__2dv1) + QLM(e12()%v1__2, dv2__12dv1) + QLM(e12()%v2__12, dv1__2dv1) + QLM(e22()%v2__2, dv2__12dv1) + QLM(e22()%v2__12, dv2__2dv1);
            arma::mat dG_22__1dv2 = dv2__12dv2 + QLM(e11()%v1__2, dv1__12dv2) + QLM(e11()%v1__12, dv1__2dv2) + QLM(e12()%v2__2, dv1__12dv2) + QLM(e12()%v1__12, dv2__2dv2) + QLM(e12()%v1__2, dv2__12dv2) + QLM(e12()%v2__12, dv1__2dv2) + QLM(e22()%v2__2, dv2__12dv2) + QLM(e22()%v2__12, dv2__2dv2);
            arma::mat dG_22__2dv1 = dv2__22dv1 + QLM(e11()%v1__2, dv1__22dv1) + QLM(e11()%v1__22, dv1__2dv1) + QLM(e12()%v2__2, dv1__22dv1) + QLM(e12()%v1__22, dv2__2dv1) + QLM(e12()%v1__2, dv2__22dv1) + QLM(e12()%v2__22, dv1__2dv1) + QLM(e22()%v2__2, dv2__22dv1) + QLM(e22()%v2__22, dv2__2dv1);
            arma::mat dG_22__2dv2 = dv2__22dv2 + QLM(e11()%v1__2, dv1__22dv2) + QLM(e11()%v1__22, dv1__2dv2) + QLM(e12()%v2__2, dv1__22dv2) + QLM(e12()%v1__22, dv2__2dv2) + QLM(e12()%v1__2, dv2__22dv2) + QLM(e12()%v2__22, dv1__2dv2) + QLM(e22()%v2__2, dv2__22dv2) + QLM(e22()%v2__22, dv2__2dv2);

            arma::mat dG_11dz  = QLM(z_1, DD1);
            arma::mat dG_22dz  = QLM(z_2, DD2);
            arma::mat dG_12dz  = QLM(z_1/2, DD2) + QLM(z_2/2, DD1);

            arma::mat dn11dv1 = D*(QLM(H1111, dG_11dv1) + QLM(2*H1112, dG_12dv1) + QLM(H1122, dG_22dv1));
            arma::mat dn11dv2 = D*(QLM(H1111, dG_11dv2) + QLM(2*H1112, dG_12dv2) + QLM(H1122, dG_22dv2));
            arma::mat dn11dz  = D*(QLM(H1111, dG_11dz)  + QLM(2*H1112, dG_12dz)  + QLM(H1122, dG_22dz));

            arma::mat dn12dv1 = D*(QLM(H1112, dG_11dv1) + QLM(2*H1221, dG_12dv1) + QLM(H1222, dG_22dv1));
            arma::mat dn12dv2 = D*(QLM(H1112, dG_11dv2) + QLM(2*H1221, dG_12dv2) + QLM(H1222, dG_22dv2));
            arma::mat dn12dz  = D*(QLM(H1112, dG_11dz)  + QLM(2*H1221, dG_12dz)  + QLM(H1222, dG_22dz));

            arma::mat dn22dv1 = D*(QLM(H1122, dG_11dv1) + QLM(2*H1222, dG_12dv1) + QLM(H2222, dG_22dv1));
            arma::mat dn22dv2 = D*(QLM(H1122, dG_11dv2) + QLM(2*H1222, dG_12dv2) + QLM(H2222, dG_22dv2));
            arma::mat dn22dz  = D*(QLM(H1122, dG_11dz)  + QLM(2*H1222, dG_12dz)  + QLM(H2222, dG_22dz));

            arma::mat dn11__1dv1 = D*(QLM(H1111, dG_11__1dv1) + QLM(2*H1112, dG_12__1dv1) + QLM(H1122, dG_22__1dv1));
            arma::mat dn11__1dv2 = D*(QLM(H1111, dG_11__1dv2) + QLM(2*H1112, dG_12__1dv2) + QLM(H1122, dG_22__1dv2));
            arma::mat dn11__1dz  = D*(QLM(H1111, dG_11__1dz)  + QLM(2*H1112, dG_12__1dz)  + QLM(H1122, dG_22__1dz));

            arma::mat dn12__1dv1 = D*(QLM(H1112, dG_11__1dv1) + QLM(2*H1221, dG_12__1dv1) + QLM(H1222, dG_22__1dv1));
            arma::mat dn12__1dv2 = D*(QLM(H1112, dG_11__1dv2) + QLM(2*H1221, dG_12__1dv2) + QLM(H1222, dG_22__1dv2));
            arma::mat dn12__1dz  = D*(QLM(H1112, dG_11__1dz)  + QLM(2*H1221, dG_12__1dz)  + QLM(H1222, dG_22__1dz));

            arma::mat dn12__2dv1 = D*(QLM(H1112, dG_11__2dv1) + QLM(2*H1221, dG_12__2dv1) + QLM(H1222, dG_22__2dv1));
            arma::mat dn12__2dv2 = D*(QLM(H1112, dG_11__2dv2) + QLM(2*H1221, dG_12__2dv2) + QLM(H1222, dG_22__2dv2));
            arma::mat dn12__2dz  = D*(QLM(H1112, dG_11__2dz)  + QLM(2*H1221, dG_12__2dz)  + QLM(H1222, dG_22__2dz));

            arma::mat dn22__2dv1 = D*(QLM(H1122, dG_11__2dv1) + QLM(2*H1222, dG_12__2dv1) + QLM(H2222, dG_22__2dv1));
            arma::mat dn22__2dv2 = D*(QLM(H1122, dG_11__2dv2) + QLM(2*H1222, dG_12__2dv2) + QLM(H2222, dG_22__2dv2));
            arma::mat dn22__2dz  = D*(QLM(H1122, dG_11__2dz)  + QLM(2*H1222, dG_12__2dz)  + QLM(H2222, dG_22__2dz));

            arma::mat dwdz    = QLM((-e11()%z_1 - e12()%z_2), DD1) + QLM((-e12()%z_1 - e22()%z_2), DD2);
            arma::mat dw_1dz  = QLM((e11()%v1__1 + e12()%v2__1 - 1), DD1) + QLM((e12()%v1__1 + e22()%v2__1)    , DD2);
            arma::mat dw_2dz  = QLM((e11()%v1__2 + e12()%v2__2)    , DD1) + QLM((e12()%v1__2 + e22()%v2__2 - 1), DD2);
            arma::mat dw_1dv1 = QLM(e11()%z_1, dv1__1dv1) + QLM(e12()%z_2, dv2__1dv1);
            arma::mat dw_2dv1 = QLM(e12()%z_1, dv1__2dv1) + QLM(e22()%z_2, dv2__2dv1);
            arma::mat dw_1dv2 = QLM(e11()%z_1, dv1__1dv2) + QLM(e12()%z_2, dv2__1dv2);
            arma::mat dw_2dv2 = QLM(e12()%z_1, dv1__2dv2) + QLM(e22()%z_2, dv2__2dv2);
            arma::mat dw1dv1  = QLM(e11(), dw_1dv1) + QLM(e12(), dw_2dv1);
            arma::mat dw2dv1  = QLM(e12(), dw_1dv1) + QLM(e22(), dw_2dv1);
            arma::mat dw1dv2  = QLM(e11(), dw_1dv2) + QLM(e12(), dw_2dv2);
            arma::mat dw2dv2  = QLM(e12(), dw_1dv2) + QLM(e22(), dw_2dv2);
            arma::mat dw1dz   = QLM(e11(), dw_1dz)  + QLM(e12(), dw_2dz);
            arma::mat dw2dz   = QLM(e12(), dw_1dz)  + QLM(e22(), dw_2dz);
            arma::mat ds1dv1  = QLM(surface1%e11(), dv1__1dv1) + QLM(surface1%e12(), dv2__1dv1) + QLM(surface2%e11(), dv1__2dv1) + QLM(surface2%e12(), dv2__2dv1) + QLM(pressure, dw1dv1);
            arma::mat ds2dv1  = QLM(surface1%e12(), dv1__1dv1) + QLM(surface1%e22(), dv2__1dv1) + QLM(surface2%e12(), dv1__2dv1) + QLM(surface2%e22(), dv2__2dv1) + QLM(pressure, dw2dv1);
            arma::mat ds1dv2  = QLM(surface1%e11(), dv1__1dv2) + QLM(surface1%e12(), dv2__1dv2) + QLM(surface2%e11(), dv1__2dv2) + QLM(surface2%e12(), dv2__2dv2) + QLM(pressure, dw1dv2);
            arma::mat ds2dv2  = QLM(surface1%e12(), dv1__1dv2) + QLM(surface1%e22(), dv2__1dv2) + QLM(surface2%e12(), dv1__2dv2) + QLM(surface2%e22(), dv2__2dv2) + QLM(pressure, dw2dv2);
            arma::mat ds1dz   = QLM(pressure, dw1dz);
            arma::mat ds2dz   = QLM(pressure, dw2dz);
            arma::mat ds3dz   = QLM(surface1, DD1) + QLM(surface2, DD2) + QLM(pressure, dwdz);

            arma::mat dL111dv1 = QLM(e11(), dv1__11dv1) + QLM(e12(), dv2__11dv1);
            arma::mat dL111dv2 = QLM(e11(), dv1__11dv2) + QLM(e12(), dv2__11dv2);
            arma::mat dL112dv1 = QLM(e11(), dv1__12dv1) + QLM(e12(), dv2__12dv1);
            arma::mat dL112dv2 = QLM(e11(), dv1__12dv2) + QLM(e12(), dv2__12dv2);
            arma::mat dL122dv1 = QLM(e11(), dv1__22dv1) + QLM(e12(), dv2__22dv1);
            arma::mat dL122dv2 = QLM(e11(), dv1__22dv2) + QLM(e12(), dv2__22dv2);
            arma::mat dL211dv1 = QLM(e12(), dv1__11dv1) + QLM(e22(), dv2__11dv1);
            arma::mat dL211dv2 = QLM(e12(), dv1__11dv2) + QLM(e22(), dv2__11dv2);
            arma::mat dL212dv1 = QLM(e12(), dv1__12dv1) + QLM(e22(), dv2__12dv1);
            arma::mat dL212dv2 = QLM(e12(), dv1__12dv2) + QLM(e22(), dv2__12dv2);
            arma::mat dL222dv1 = QLM(e12(), dv1__22dv1) + QLM(e22(), dv2__22dv1);
            arma::mat dL222dv2 = QLM(e12(), dv1__22dv2) + QLM(e22(), dv2__22dv2);

            arma::mat dsaedz  = QLM(e11(), dG_11dz)  + QLM(2*e12(), dG_12dz)  + QLM(e22(), dG_22dz);
            arma::mat dsaedv1 = QLM(e11(), dG_11dv1) + QLM(2*e12(), dG_12dv1) + QLM(e22(), dG_22dv1);
            arma::mat dsaedv2 = QLM(e11(), dG_11dv2) + QLM(2*e12(), dG_12dv2) + QLM(e22(), dG_22dv2);

            arma::mat A(nxy3, nxy3, arma::fill::zeros);

            A.submat(   0,    0,  nxy-1,  nxy-1) = (dn11__1dv1 + dn12__2dv1 + QLM(L111, dn11dv1) + QLM(2*L112, dn12dv1) + QLM(L122, dn22dv1) + QLM(sae, ds1dv1) + QLM(s1, dsaedv1) + QLM(n11, dL111dv1) + QLM(2*n12, dL112dv1) + QLM(n22, dL122dv1))/D;
            A.submat(   0,  nxy,  nxy-1, nxy2-1) = (dn11__1dv2 + dn12__2dv2 + QLM(L111, dn11dv2) + QLM(2*L112, dn12dv2) + QLM(L122, dn22dv2) + QLM(sae, ds1dv2) + QLM(s1, dsaedv2) + QLM(n11, dL111dv2) + QLM(2*n12, dL112dv2) + QLM(n22, dL122dv2))/D;
            A.submat(   0, nxy2,  nxy-1, nxy3-1) = (dn11__1dz  + dn12__2dz  + QLM(L111, dn11dz)  + QLM(2*L112, dn12dz)  + QLM(L122, dn22dz)  + QLM(sae, ds1dz)  + QLM(s1, dsaedz))/D;

            A.submat( nxy,    0, nxy2-1,  nxy-1) = (dn12__1dv1 + dn22__2dv1 + QLM(L211, dn11dv1) + QLM(2*L212, dn12dv1) + QLM(L222, dn22dv1) + QLM(sae, ds2dv1) + QLM(s2, dsaedv1) + QLM(n11, dL211dv1) + QLM(2*n12, dL212dv1) + QLM(n22, dL222dv1))/D;
            A.submat( nxy,  nxy, nxy2-1, nxy2-1) = (dn12__1dv2 + dn22__2dv2 + QLM(L211, dn11dv2) + QLM(2*L212, dn12dv2) + QLM(L222, dn22dv2) + QLM(sae, ds2dv2) + QLM(s2, dsaedv2) + QLM(n11, dL211dv2) + QLM(2*n12, dL212dv2) + QLM(n22, dL222dv2))/D;
            A.submat( nxy, nxy2, nxy2-1, nxy3-1) = (dn12__1dz  + dn22__2dz  + QLM(L211, dn11dz)  + QLM(2*L212, dn12dz)  + QLM(L222, dn22dz)  + QLM(sae, ds2dz)  + QLM(s2, dsaedz))/D;

            A.submat(nxy2,    0, nxy3-1,  nxy-1) = (QLM(z__11, dn11dv1) + QLM(2*z__12, dn12dv1) + QLM(z__22, dn22dv1) + QLM(s3, dsaedv1))/D;
            A.submat(nxy2,  nxy, nxy3-1, nxy2-1) = (QLM(z__11, dn11dv2) + QLM(2*z__12, dn12dv2) + QLM(z__22, dn22dv2) + QLM(s3, dsaedv2))/D;
            A.submat(nxy2, nxy2, nxy3-1, nxy3-1) = (QLM(z__11, dn11dz)  + QLM(2*z__12, dn12dz)  + QLM(z__22, dn22dz)  + QLM(s3, dsaedz) + QLM(sae, ds3dz) + QLM(n11, dz__11dz) + QLM(2*n12, dz__12dz) + QLM(n22, dz__22dz))/D;

            for (size_t j = 1; j < ny-1; j++)
            {
                // BC west
                v1BoundaryWestEastNonlinear(v1.westBC, v1.west(j), 0, j, A, v1.r1West, v1.r2West, h_1s1_west(j), h_1s2_west(j));
                v2BoundaryWestEastNonlinear(v2.westBC, v2.west(j), 0, j, A, v2.r1West, v2.r2West, h_1s1_west(j), h_1s2_west(j), h_2s2_west(j));
                n11BoundaryNonlinear(n11.westBC, n11.west(j), 0, j, A, h_11s_west(j), z_1(0, j), z_2(0, j));
                n12BoundaryNonlinear(n12.westBC, n12.west(j), 0, j, A, h_11s_west(j), h_12s_west(j), h_22s_west(j), z_1(0, j), z_2(0, j));
                zBoundaryNonlinear(z.westBC, z.west(j), 0, j, A, z.r1West, z.r2West, z_1(0, j), z_2(0, j), h_1s1_west(j), h_1s2_west(j));

                // BC east
                size_t i = nx-1;
                v1BoundaryWestEastNonlinear(v1.eastBC, v1.east(j), i, j, A, v1.r1East, v1.r2East, h_1s1_east(j), h_1s2_east(j));
                v2BoundaryWestEastNonlinear(v2.eastBC, v2.east(j), i, j, A, v2.r1East, v2.r2East, h_1s1_east(j), h_1s2_east(j), h_2s2_east(j));
                n11BoundaryNonlinear(n11.eastBC, n11.east(j), i, j, A, h_11s_east(j), z_1(i, j), z_2(i, j));
                n12BoundaryNonlinear(n12.eastBC, n12.east(j), i, j, A, h_11s_east(j), h_12s_east(j), h_22s_east(j), z_1(i, j), z_2(i, j));
                zBoundaryNonlinear(z.eastBC, z.east(j), i, j, A, z.r1East, z.r2East, z_1(i, j), z_2(i, j), h_1s1_east(j), h_1s2_east(j));
            }
            for (size_t i = 1; i < nx-1; i++)
            {
                // BC south
                v1BoundarySouthNorthNonlinear(v1.southBC, v1.south(i), i, 0, A, v1.r1South, v1.r2South, h_1s1_south(i), h_2s1_south(i), h_2s2_south(i));
                v2BoundarySouthNorthNonlinear(v2.southBC, v2.south(i), i, 0, A, v2.r1South, v2.r2South, h_2s1_south(i), h_2s2_south(i));
                n21BoundaryNonlinear(n12.southBC, n12.south(i), i, 0, A, h_11s_south(i), h_21s_south(i), h_22s_south(i), z_1(i, 0), z_2(i, 0));
                n22BoundaryNonlinear(n22.southBC, n22.south(i), i, 0, A, h_22s_south(i), z_1(i, 0), z_2(i, 0));
                zBoundaryNonlinear(z.southBC, z.south(i), i, 0, A, z.r1South, z.r2South, z_1(i, 0), z_2(i, 0), h_2s1_south(i), h_2s2_south(i));

                // BC north
                size_t j = ny-1;
                v1BoundarySouthNorthNonlinear(v1.northBC, v1.north(i), i, j, A, v1.r1North, v1.r2North, h_1s1_north(i), h_2s1_north(i), h_2s2_north(i));
                v2BoundarySouthNorthNonlinear(v2.northBC, v2.north(i), i, j, A, v2.r1North, v2.r2North, h_2s1_north(i), h_2s2_north(i));
                n21BoundaryNonlinear(n12.northBC, n12.north(i), i, j, A, h_11s_north(i), h_21s_north(i), h_22s_north(i), z_1(i, j), z_2(i, j));
                n22BoundaryNonlinear(n22.northBC, n22.north(i), i, j, A, h_22s_north(i), z_1(i, j), z_2(i, j));
                zBoundaryNonlinear(z.northBC, z.north(i), i, j, A, z.r1North, z.r2North, z_1(i, j), z_2(i, j), h_2s1_north(i), h_2s2_north(i));
            }
            
            // BC south-west corner (i = 0, j = 0)
            v1BoundarySouthNorthNonlinear(v1.southBC, v1.south(0), 0, 0, A, v1.r1South, v1.r2South, h_1s1_south(0), h_2s1_south(0), h_2s2_south(0));
            v2BoundarySouthNorthNonlinear(v2.southBC, v2.south(0), 0, 0, A, v2.r1South, v2.r2South, h_2s1_south(0), h_2s2_south(0));
            n21BoundaryNonlinear(n12.southBC, n12.south(0), 0, 0, A, h_11s_south(0), h_21s_south(0), h_22s_south(0), z_1(0, 0), z_2(0, 0));
            n22BoundaryNonlinear(n22.southBC, n22.south(0), 0, 0, A, h_22s_south(0), z_1(0, 0), z_2(0, 0));
            zBoundaryNonlinear(z.southBC, z.south(0), 0, 0, A, z.r1South, z.r2South, z_1(0, 0), z_2(0, 0), h_2s1_south(0), h_2s2_south(0));

            // BC north-west corner (i = 0, j = ny-1)
            size_t j = ny-1;
            v1BoundaryWestEastNonlinear(v1.westBC, v1.west(j), 0, j, A, v1.r1West, v1.r2West, h_1s1_west(j), h_1s2_west(j));
            v2BoundaryWestEastNonlinear(v2.westBC, v2.west(j), 0, j, A, v2.r1West, v2.r2West, h_1s1_west(j), h_1s2_west(j), h_2s2_west(j));
            n11BoundaryNonlinear(n11.westBC, n11.west(j), 0, j, A, h_11s_west(j), z_1(0, j), z_2(0, j));
            n12BoundaryNonlinear(n12.westBC, n12.west(j), 0, j, A, h_11s_west(j), h_12s_west(j), h_22s_west(j), z_1(0, j), z_2(0, j));
            zBoundaryNonlinear(z.westBC, z.west(j), 0, j, A, z.r1West, z.r2West, z_1(0, j), z_2(0, j), h_1s1_west(j), h_1s2_west(j));

            // BC south-east corner (i = nx-1, j = 0)
            size_t i = nx-1;
            v1BoundaryWestEastNonlinear(v1.eastBC, v1.east(0), i, 0, A, v1.r1East, v1.r2East, h_1s1_east(0), h_1s2_east(0));
            v2BoundaryWestEastNonlinear(v2.eastBC, v2.east(0), i, 0, A, v2.r1East, v2.r2East, h_1s1_east(0), h_1s2_east(0), h_2s2_east(0));
            n11BoundaryNonlinear(n11.eastBC, n11.east(0), i, 0, A, h_11s_east(0), z_1(i, 0), z_2(i, 0));
            n12BoundaryNonlinear(n12.eastBC, n12.east(0), i, 0, A, h_11s_east(0), h_12s_east(0), h_22s_east(0), z_1(i, 0), z_2(i, 0));
            zBoundaryNonlinear(z.eastBC, z.east(0), i, 0, A, z.r1East, z.r2East, z_1(i, 0), z_2(i, 0), h_1s1_east(0), h_1s2_east(0));

            // BC north-east corner (i = nx-1, j = ny-1)
            i = nx-1;
            j = ny-1;
            v1BoundarySouthNorthNonlinear(v1.northBC, v1.north(i), i, j, A, v1.r1North, v1.r2North, h_1s1_north(i), h_2s1_north(i), h_2s2_north(i));
            v2BoundarySouthNorthNonlinear(v2.northBC, v2.north(i), i, j, A, v2.r1North, v2.r2North, h_2s1_north(i), h_2s2_north(i));
            n21BoundaryNonlinear(n12.northBC, n12.north(i), i, j, A, h_11s_north(i), h_21s_north(i), h_22s_north(i), z_1(i, j), z_2(i, j));
            n22BoundaryNonlinear(n22.northBC, n22.north(i), i, j, A, h_22s_north(i), z_1(i, j), z_2(i, j));
            zBoundaryNonlinear(z.northBC, z.north(i), i, j, A, z.r1North, z.r2North, z_1(i, j), z_2(i, j), h_2s1_north(i), h_2s2_north(i));
            
            arma::vec res1 = b(arma::span(   0,  nxy-1));
            arma::vec res2 = b(arma::span( nxy, nxy2-1));
            arma::vec res3 = b(arma::span(nxy2, nxy3-1));

            double res1MAX = max(arma::abs(res1))/pMAX;
            double res2MAX = max(arma::abs(res2))/pMAX;
            double res3MAX = max(arma::abs(res3))/pMAX;
            double res1RMS = sqrt(sum(res1%res1)/nxy)/pRMS;
            double res2RMS = sqrt(sum(res2%res2)/nxy)/pRMS;
            double res3RMS = sqrt(sum(res3%res3)/nxy)/pRMS;
            std::println("Res 1\t\t2\t\t3");
            std::println("Max {0:4.2e}\t{1:4.2e}\t{2:4.2e}", res1MAX, res2MAX, res3MAX);
            std::println("RMS {0:4.2e}\t{1:4.2e}\t{2:4.2e}", res1RMS, res2RMS, res3RMS);
            if (res1RMS < residualTarget && res2RMS < residualTarget && res3RMS < residualTarget)
                break;

            arma::vec dv = solve(A, b);

            double omega = armijoNonlinear(dv, surface1, surface2, pressure);
            std::println("omega = {0:4.2f}", omega);
            v1.arma::mat::operator-=(reshape(omega*dv(arma::span(   0,  nxy-1)), nx, ny));
            v2.arma::mat::operator-=(reshape(omega*dv(arma::span( nxy, nxy2-1)), nx, ny));
            z .arma::mat::operator-=(reshape(omega*dv(arma::span(nxy2, nxy3-1)), nx, ny));
            
            v1_1  = D1*v1;
            v1_2  = v1*D2_t;
            v2_1  = D1*v2;
            v2_2  = v2*D2_t;
            v1__1 = v1_1 - g111()%v1 - g211()%v2;
            v1__2 = v1_2 - g112()%v1 - g212()%v2;
            v2__1 = v2_1 - g112()%v1 - g212()%v2;
            v2__2 = v2_2 - g122()%v1 - g222()%v2;
            z_1 = D1*z;
            z_2 = z*D2_t;

            gamma_11 =  v1__1         + (z_1%z_1 + e11()%v1__1%v1__1 + 2*e12()%v1__1%v2__1                 + e22()%v2__1%v2__1)/2;
            gamma_12 = (v1__2 + v2__1 +  z_1%z_2 + e11()%v1__1%v1__2 +   e12()%(v1__1%v2__2 + v2__1%v1__2) + e22()%v2__1%v2__2)/2;
            gamma_22 =  v2__2         + (z_2%z_2 + e11()%v1__2%v1__2 + 2*e12()%v1__2%v2__2                 + e22()%v2__2%v2__2)/2;

            n11 = D*(H1111%gamma_11 + 2*H1112%gamma_12 + H1122%gamma_22);
            n12 = D*(H1112%gamma_11 + 2*H1221%gamma_12 + H1222%gamma_22);
            n22 = D*(H1122%gamma_11 + 2*H1222%gamma_12 + H2222%gamma_22);
        }
    }
}

double Membrane::armijoNonlinear(arma::vec dv, arma::mat &surface1, arma::mat &surface2, arma::mat &pressure)
{
    double omega = 1;
    arma::vec v = join_vert(vectorise(arma::mat(v1)), vectorise(arma::mat(v2)), vectorise(arma::mat(z)));
    double T = residualLevelFunctionNonlinear(v, surface1, surface2, pressure);
    double T_dx_min = 1e10;
    for (size_t m = 0; m < 20; m++)
    {
        double alpha = pow(0.5, m);
        double T_dx = residualLevelFunctionNonlinear(v - alpha*dv, surface1, surface2, pressure);
        if (T_dx <= (1 - alpha/2)*T)
            if (T_dx < T_dx_min)
            {
                omega = alpha;
                T_dx_min = T_dx;
            }
    }
    return omega;
}

double Membrane::residualLevelFunctionNonlinear(arma::vec V, arma::mat &surface1, arma::mat &surface2, arma::mat &pressure)
{
    arma::mat V1 = reshape(V(arma::span(    0,   nxy-1)), nx, ny);
    arma::mat V2 = reshape(V(arma::span(  nxy, 2*nxy-1)), nx, ny);
    arma::mat Z  = reshape(V(arma::span(2*nxy, 3*nxy-1)), nx, ny);
    arma::mat V1_1 = D1*V1;
    arma::mat V1_2 = V1*D2.t();
    arma::mat V2_1 = D1*V2;
    arma::mat V2_2 = V2*D2.t();
    arma::mat V1__1 = V1_1 - g111()%V1 - g211()%V2;
    arma::mat V1__2 = V1_2 - g112()%V1 - g212()%V2;
    arma::mat V2__1 = V2_1 - g112()%V1 - g212()%V2;
    arma::mat V2__2 = V2_2 - g122()%V1 - g222()%V2;
    arma::mat Z_1  = D1*Z;
    arma::mat Z_2  = Z*D2.t();

    arma::mat V1_11 = D11*V1;
    arma::mat V1_12 = D1*V1*D2.t();
    arma::mat V1_22 = V1*D22.t();
    arma::mat V2_11 = D11*V2;
    arma::mat V2_12 = D1*V2*D2.t();
    arma::mat V2_22 = V2*D22.t();
    arma::mat Z_11 = D11*Z;
    arma::mat Z_12 = D1*Z*D2.t();
    arma::mat Z_22 = Z*D22.t();

    arma::mat V1__11 = V1_11 - 3*g111()%V1_1 - g211()%(V1_2 + 2*V2_1) + Psi1111%V1 + Psi1112%V2;
    arma::mat V2__22 = V2_22 - 3*g222()%V2_2 - g122()%(V2_1 + 2*V1_2) + Psi2221%V1 + Psi2222%V2;
    arma::mat V1__12 = V1_12 - 2*g112()%V1_1 - (g111() + g212())%V1_2 - g212()%V2_1 - g211()%V2_2 + Psi1121%V1 + Psi1122%V2;
    arma::mat V2__12 = V2_12 - 2*g212()%V2_2 - (g112() + g222())%V2_1 - g122()%V1_1 - g112()%V1_2 + Psi2121%V1 + Psi2122%V2;
    arma::mat V1__22 = V1_22 -   g122()%V1_1 - (2*g112() + g222())%V1_2 - 2*g212()%V2_2 + Psi1221%V1 + Psi1222%V2;
    arma::mat V2__11 = V2_11 -   g211()%V2_2 - (2*g212() + g111())%V2_1 - 2*g112()%V1_1 + Psi2112%V1 + Psi2111%V2;
    arma::mat Z__11  = Z_11  -   g111()%Z_1  - g211()%Z_2;
    arma::mat Z__12  = Z_12  -   g112()%Z_1  - g212()%Z_2;
    arma::mat Z__22  = Z_22  -   g122()%Z_1  - g222()%Z_2;

    arma::mat GAMMA_11 =  V1__1         + (Z_1%Z_1 + e11()%V1__1%V1__1 + 2*e12()%V1__1%V2__1                 + e22()%V2__1%V2__1)/2;
    arma::mat GAMMA_12 = (V1__2 + V2__1 +  Z_1%Z_2 + e11()%V1__1%V1__2 +   e12()%(V1__1%V2__2 + V2__1%V1__2) + e22()%V2__1%V2__2)/2;
    arma::mat GAMMA_22 =  V2__2         + (Z_2%Z_2 + e11()%V1__2%V1__2 + 2*e12()%V1__2%V2__2                 + e22()%V2__2%V2__2)/2;

    arma::mat Z_1__12 = Z_1%Z__12;
    arma::mat Z_2__12 = Z_2%Z__12;
    arma::mat Z_112 = Z_1__12 + Z_2%Z__11;
    arma::mat Z_122 = Z_1%Z__22 + Z_2__12;

    arma::mat G_11__1 = V1__11 + Z_1%Z__11         + e11()%V1__1%V1__11                  + e12()%(V1__1%V2__11 + V2__1%V1__11)                               + e22()%V2__1%V2__11;
    arma::mat G_11__2 = V1__12 + Z_1__12           + e11()%V1__1%V1__12                  + e12()%(V1__1%V2__12 + V2__1%V1__12)                               + e22()%V2__1%V2__12;
    arma::mat G_12__1 = (V1__12 + V2__11   + Z_112 + e11()%(V1__1%V1__12 + V1__2%V1__11) + e12()%(V2__1%V1__12 + V1__2%V2__11 + V1__1%V2__12 + V2__2%V1__11) + e22()%(V2__1%V2__12 + V2__2%V2__11))/2;
    arma::mat G_12__2 = (V1__22 + V2__12   + Z_122 + e11()%(V1__1%V1__22 + V1__2%V1__12) + e12()%(V2__1%V1__22 + V1__2%V2__12 + V1__1%V2__22 + V2__2%V1__12) + e22()%(V2__1%V2__22 + V2__2%V2__12))/2;
    arma::mat G_22__1 = V2__12 + Z_2__12           + e11()%V1__2%V1__12                  + e12()%(V1__2%V2__12 + V2__2%V1__12)                               + e22()%V2__2%V2__12;
    arma::mat G_22__2 = V2__22 + Z_2%Z__22         + e11()%V1__2%V1__22                  + e12()%(V1__2%V2__22 + V2__2%V1__22)                               + e22()%V2__2%V2__22;

    arma::mat N11 = D*(H1111%GAMMA_11 + 2*H1112%GAMMA_12 + H1122%GAMMA_22);
    arma::mat N12 = D*(H1112%GAMMA_11 + 2*H1221%GAMMA_12 + H1222%GAMMA_22);
    arma::mat N22 = D*(H1122%GAMMA_11 + 2*H1222%GAMMA_12 + H2222%GAMMA_22);

    arma::mat N11__1 = D*(H1111%G_11__1 + 2*H1112%G_12__1 + H1122%G_22__1);
    arma::mat N12__1 = D*(H1112%G_11__1 + 2*H1221%G_12__1 + H1222%G_22__1);
    arma::mat N12__2 = D*(H1112%G_11__2 + 2*H1221%G_12__2 + H1222%G_22__2);
    arma::mat N22__2 = D*(H1122%G_11__2 + 2*H1222%G_12__2 + H2222%G_22__2);

    arma::mat L111 = e11()%V1__11 + e12()%V2__11;
    arma::mat L112 = e11()%V1__12 + e12()%V2__12;
    arma::mat L122 = e11()%V1__22 + e12()%V2__22;
    arma::mat L211 = e12()%V1__11 + e22()%V2__11;
    arma::mat L212 = e12()%V1__12 + e22()%V2__12;
    arma::mat L222 = e12()%V1__22 + e22()%V2__22;

    arma::mat SAE = 1 + e11()%GAMMA_11 + 2*e12()%GAMMA_12 + e22()%GAMMA_22;

    arma::mat W   =-(e_11()%Z_1%Z_1 + 2*e_12()%Z_1%Z_2 + e_22()%Z_2%Z_2)/2;
    arma::mat W_1 = Z_1%(e_11()%V1__1 + e_12()%V2__1 - 1) + Z_2%(e_12()%V1__1 + e_22()%V2__1);
    arma::mat W_2 = Z_1%(e_11()%V1__2 + e_12()%V2__2) + Z_2%(e_12()%V1__2 + e_22()%V2__2 - 1);

    arma::mat W1  = e_11()%W_1 + e_12()%W_2;
    arma::mat W2  = e_12()%W_1 + e_22()%W_2;
    
    arma::mat S1 = surface1%(e_11()%V1__1 + e_12()%V2__1 + 1) + surface2%(e_11()%V1__2 + e_12()%V2__2)     + pressure%W1;
    arma::mat S2 = surface1%(e_12()%V1__1 + e_22()%V2__1)     + surface2%(e_12()%V1__2 + e_22()%V2__2 + 1) + pressure%W2;
    arma::mat S3 = surface1%Z_1 + surface2%Z_2 + pressure%(1 + W);

    arma::vec f = join_vert(vectorise(N11__1 + N12__2 + L111%N11 + 2*L112%N12 + L122%N22 + SAE%S1)/D,
                            vectorise(N12__1 + N22__2 + L211%N11 + 2*L212%N12 + L222%N22 + SAE%S2)/D,
                            vectorise(Z__11%N11 + 2*Z__12%N12 + Z__22%N22 + SAE%S3)/D);
    return dot(f, f);
}