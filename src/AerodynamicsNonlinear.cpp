#include "Aerodynamics.hpp"

void Aerodynamics::nonlinear()
{
    analysis = Analysis::nonlinear;
    // Influence of the wing surfaces on each other
    bw.set_size(wings.size(), wings.size());
    for (size_t sD = 0; sD < wings.size(); sD++)
    {
        wings[sD]->analysis = Analysis::nonlinear;
        for (size_t tD = 0; tD < wings.size(); tD++)
            if (sD != tD)
            {
                std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                {
                    gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                    x1_gl(ii) =-gl_x[ii].x();
                }
                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                {
                    gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                    x2_gl(jj) =-gl_y[jj].x();
                }
                arma::mat Tx = Lagrange::interpolationMatrix(wings[sD]->x1, x1_gl);
                arma::mat Ty = Lagrange::interpolationMatrix(wings[sD]->x2, x2_gl);
                arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, wings[sD]->z, x1_gl, x2_gl);
                arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
                arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
                arma::mat dz_gldx1 = D1_gl*z_gl;
                arma::mat dz_gldx2 = z_gl*D2_gl.t();
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
                arma::cube e_c_gl = MetricCo(J_gl);
                arma::cube ec_gl  = MetricContra(e_c_gl);
                arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
                arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
                arma::mat xC = wings[tD]->xC;
                arma::mat yC = wings[tD]->yC;
                arma::mat zC = wings[tD]->zC;
                bw(tD, sD).set_size(wings[tD]->nxy, wings[sD]->nxy, wings[sD]->con);
                for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                    for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                    {
                        size_t k = i+j*wings[tD]->nx;
                        for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                            for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                            {
                                arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
                                double r3 = pow(norm(r), 3);
                                arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                                                dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                                                dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
                                arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                                                {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                                                {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
                                for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                {
                                    double  t2 = boost::math::chebyshev_t(q, x2_gl(jj));
                                    double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                    for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                    {
                                        double  t1 = boost::math::chebyshev_t(p, x1_gl(ii));
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                        arma::vec::fixed<2> dmudxi = {dt1 *  t2, t1 * dt2};
                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                        arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                                        arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3;
                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k));
                                    }
                                }
                            } 
                    }
                for (const Wake* w:wings[sD]->wakes)
                    for (size_t c = 0; c < 4; c++)
                        if (wings[sD]->chi[c] == w->chi)
                        {
                            std::vector<fastgl::QuadPair> gl_x_w(wings[sD]->nx), gl_y_w(wings[sD]->ny);
                            arma::vec x1_gl_w(wings[sD]->nx), x2_gl_w(wings[sD]->ny);
                            for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                            {
                                gl_x_w[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                                x1_gl_w(ii) =-gl_x_w[ii].x();
                            }
                            for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                            {
                                gl_y_w[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                                x2_gl_w(jj) =-gl_y_w[jj].x();
                            }
                            arma::vec x1_gauss = Chebyshev::gaussLobatto(wings[sD]->nx);
                            arma::vec x2_gauss = Chebyshev::gaussLobatto(wings[sD]->ny);
                            arma::mat D1_gauss = Chebyshev::derivativeMatrix(x1_gauss, Derivative::first);
                            arma::mat D2_gauss = Chebyshev::derivativeMatrix(x2_gauss, Derivative::first);
                            arma::mat dxdx1_gauss = D1_gauss * wings[sD]->x;
                            arma::mat dydx1_gauss = D1_gauss * wings[sD]->y;
                            arma::mat dzdx1_gauss = D1_gauss * wings[sD]->z;
                            arma::mat dxdx2_gauss = wings[sD]->x * D2_gauss.t();
                            arma::mat dydx2_gauss = wings[sD]->y * D2_gauss.t();
                            arma::mat dzdx2_gauss = wings[sD]->z * D2_gauss.t();
                            arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                            arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                            if (c == 0) // South
                            {
                                arma::vec xw    = T1_gauss * wings[sD]->x.col(0);
                                arma::vec yw    = T1_gauss * wings[sD]->y.col(0);
                                arma::vec zw    = T1_gauss * wings[sD]->z.col(0);
                                arma::vec dxdx1 = T1_gauss *  dxdx1_gauss.col(0);
                                arma::vec dydx1 = T1_gauss *  dydx1_gauss.col(0);
                                arma::vec dzdx1 = T1_gauss *  dzdx1_gauss.col(0);
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                {
                                    arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                            double zW = zw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))*tan(wings[sD]->alpha(0))/2;
                                            double dxdx2 =-pow(1 + x2_gl_w(jj),-2);
                                            // dydx2 = 0;
                                            double dzdx2 =-tan(wings[sD]->alpha(0))/pow(1 + x2_gl_w(jj), 2);
                                            arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t q = 0; q < wings[sD]->ny; q++)
                                            {
                                                double t2 = pow(-1, q);
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = dt1 * t2;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                       dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                      -dydx1(ii)*dxdx2})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {-n_W(1)*dzdx2,
                                                                                     n_W(0)*dzdx2 - dxdx2*n_W(2),
                                                                                     dxdx2*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k)); 
                                                    }
                                            }
                                        }
                                    }
                                }
                            }
                            else if (c == 1) // East
                            {
                                arma::vec xw    = T2_gauss * wings[sD]->x.row(wings[sD]->nx-1).t();
                                arma::vec yw    = T2_gauss * wings[sD]->y.row(wings[sD]->nx-1).t();
                                arma::vec zw    = T2_gauss * wings[sD]->z.row(wings[sD]->nx-1).t();
                                arma::vec dxdx2 = T2_gauss *  dxdx2_gauss.row(wings[sD]->nx-1).t();
                                arma::vec dydx2 = T2_gauss *  dydx2_gauss.row(wings[sD]->nx-1).t();
                                arma::vec dzdx2 = T2_gauss *  dzdx2_gauss.row(wings[sD]->nx-1).t();
                                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                {
                                    arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                            double zW = zw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))*tan(wings[sD]->alpha(0))/2;
                                            double dxdx1 = 1/pow(1 - x1_gl_w(ii), 2);
                                            // dydx1 = 0
                                            double dzdx1 = tan(wings[sD]->alpha(0))/pow(1 - x1_gl_w(ii), 2);
                                            arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                {
                                                    size_t k = i + j*wings[tD]->nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = dt2;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                    dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                    dxdx1*dydx2(jj)})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                 dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                 -dxdx1*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k));
                                                }
                                        }
                                    }
                                }
                            }
                            else if (c == 2) // North
                            {
                                arma::vec xw    = T1_gauss * wings[sD]->x.col(wings[sD]->ny-1);
                                arma::vec yw    = T1_gauss * wings[sD]->y.col(wings[sD]->ny-1);
                                arma::vec zw    = T1_gauss * wings[sD]->z.col(wings[sD]->ny-1);
                                arma::vec dxdx1 = T1_gauss *  dxdx1_gauss.col(wings[sD]->ny-1);
                                arma::vec dydx1 = T1_gauss *  dydx1_gauss.col(wings[sD]->ny-1);
                                arma::vec dzdx1 = T1_gauss *  dzdx1_gauss.col(wings[sD]->ny-1);
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                {
                                    arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                            double zW = zw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))*tan(wings[sD]->alpha(0))/2;
                                            double dxdx2 = 1/pow(1 - x2_gl_w(jj), 2);
                                            // dydx2 = 0
                                            double dzdx2 = tan(wings[sD]->alpha(0))/pow(1 - x2_gl_w(jj), 2);
                                            arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                {
                                                    size_t k = i + j*wings[tD]->nx;
                                                    arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double dmudxi = dt1;
                                                    double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                    arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                   dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                  -dydx1(ii)*dxdx2})/sqrt_A;
                                                    arma::vec::fixed<3> J_red = {-dzdx2*n_W(1),
                                                                                  dzdx2*n_W(0) - dxdx2*n_W(2),
                                                                                  dxdx2*n_W(1)};
                                                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                    arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                    arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k));
                                                }
                                        }
                                    }
                                }
                            }
                            else // West
                            {
                                arma::vec xw    = T2_gauss * wings[sD]->x.row(0).t();
                                arma::vec yw    = T2_gauss * wings[sD]->y.row(0).t();
                                arma::vec zw    = T2_gauss * wings[sD]->z.row(0).t();
                                arma::vec dxdx2 = T2_gauss *  dxdx2_gauss.row(0).t();
                                arma::vec dydx2 = T2_gauss *  dydx2_gauss.row(0).t();
                                arma::vec dzdx2 = T2_gauss *  dzdx2_gauss.row(0).t();
                                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                {
                                    arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                            double zW = zw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))*tan(wings[sD]->alpha(0))/2;
                                            double dxdx1 =-1/pow(1 + x1_gl_w(ii), 2);
                                            // dydx1 = 0
                                            double dzdx1 =-tan(wings[sD]->alpha(0))/pow(1 + x1_gl_w(ii), 2);
                                            arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                            double e_11 = dot(e_1, e_1);
                                            double e_12 = dot(e_1, e_2);
                                            double e_22 = dot(e_2, e_2);
                                            double e = e_11*e_22 - pow(e_12, 2);
                                            double e11 = e_22/e;
                                            double e12 =-e_12/e;
                                            double e22 = e_11/e;
                                            for (size_t p = 0; p < wings[sD]->nx; p++)
                                            {
                                                double t1 = pow(-1, p);
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = t1*dt2;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                        dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                        dxdx1*dydx2(jj)})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                     dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                    -dxdx1*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx, 0) -= dot(q_mu, wings[tD]->nC.row(k)); 
                                                    }
                                            }
                                        }
                                    }
                                }
                            }
                        }
            }
    }
    if (sym == Symmetry::y)
    {
        for (size_t sD = 0; sD < wings.size(); sD++)
            for (size_t tD = 0; tD < wings.size(); tD++)
                if (sD != tD)
                {
                    std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                    arma::vec x1_gl(wings[sD]->nx), x2_gl(wings[sD]->ny);
                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                    {
                        gl_x[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                        x1_gl(ii) =-gl_x[ii].x();
                    }
                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                    {
                        gl_y[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                        x2_gl(jj) =-gl_y[jj].x();
                    }
                    arma::mat Tx = Lagrange::interpolationMatrix(wings[sD]->x1, x1_gl);
                    arma::mat Ty = Lagrange::interpolationMatrix(wings[sD]->x2, x2_gl);
                    arma::mat z_gl = Lagrange::interpolation2D(Tx, Ty, wings[sD]->z, x1_gl, x2_gl);
                    arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
                    arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
                    arma::mat dz_gldx1 = D1_gl*z_gl;
                    arma::mat dz_gldx2 = z_gl*D2_gl.t();
                    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                    y_gl     =-y_gl;
                    dy_gldx1 =-dy_gldx1;
                    dy_gldx2 =-dy_gldx2;
                    arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
                    arma::cube e_c_gl = MetricCo(J_gl);
                    arma::cube ec_gl  = MetricContra(e_c_gl);
                    arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
                    arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
                    arma::mat xC = wings[tD]->xC;
                    arma::mat yC = wings[tD]->yC;
                    arma::mat zC = wings[tD]->zC;
                    for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                        for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                        {
                            size_t k = i+j*wings[tD]->nx;
                            for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                                {
                                    arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
                                    double r3 = pow(norm(r), 3);
                                    arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                                                    dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                                                    dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
                                    arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),-(dy_gldx1(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx1(ii, jj))},
                                                                    {-(dx_gldx2(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx2(ii, jj)),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                                                    {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),-(dx_gldx1(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx1(ii, jj))}};
                                    for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction
                                    {
                                        double  t2 = boost::math::chebyshev_t(q, x2_gl(jj));
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl(jj));
                                        for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction
                                        {
                                            double  t1 = boost::math::chebyshev_t(p, x1_gl(ii));
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl(ii));
                                            arma::vec::fixed<2> dmudxi = {dt1 *  t2, t1 * dt2};
                                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                            arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                                            arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3;
                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) += dot(q_mu, wings[tD]->nC.row(k));
                                        }
                                    }
                                } 
                        }
                    for (const Wake* w:wings[sD]->wakes)
                        for (size_t c = 0; c < 4; c++)
                            if (wings[sD]->chi[c] == w->chi)
                            {
                                std::vector<fastgl::QuadPair> gl_x_w(wings[sD]->nx), gl_y_w(wings[sD]->ny);
                                arma::vec x1_gl_w(wings[sD]->nx), x2_gl_w(wings[sD]->ny);
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                {
                                    gl_x_w[ii] = fastgl::GLPair(wings[sD]->nx, ii+1);
                                    x1_gl_w(ii) =-gl_x_w[ii].x();
                                }
                                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                {
                                    gl_y_w[jj] = fastgl::GLPair(wings[sD]->ny, jj+1);
                                    x2_gl_w(jj) =-gl_y_w[jj].x();
                                }
                                arma::vec x1_gauss = Chebyshev::gaussLobatto(wings[sD]->nx);
                                arma::vec x2_gauss = Chebyshev::gaussLobatto(wings[sD]->ny);
                                arma::mat D1_gauss = Chebyshev::derivativeMatrix(x1_gauss, Derivative::first);
                                arma::mat D2_gauss = Chebyshev::derivativeMatrix(x2_gauss, Derivative::first);
                                arma::mat dxdx1_gauss = D1_gauss * wings[sD]->x;
                                arma::mat dydx1_gauss = D1_gauss * wings[sD]->y;
                                arma::mat dzdx1_gauss = D1_gauss * wings[sD]->z;
                                arma::mat dxdx2_gauss = wings[sD]->x * D2_gauss.t();
                                arma::mat dydx2_gauss = wings[sD]->y * D2_gauss.t();
                                arma::mat dzdx2_gauss = wings[sD]->z * D2_gauss.t();
                                arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                                arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                                if (c == 0) // South
                                {
                                    arma::vec xw    = T1_gauss * wings[sD]->x.col(0);
                                    arma::vec yw    =-T1_gauss * wings[sD]->y.col(0);
                                    arma::vec zw    = T1_gauss * wings[sD]->z.col(0);
                                    arma::vec dxdx1 = T1_gauss *  dxdx1_gauss.col(0);
                                    arma::vec dydx1 =-T1_gauss *  dydx1_gauss.col(0);
                                    arma::vec dzdx1 = T1_gauss *  dzdx1_gauss.col(0);
                                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                    {
                                        arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                        {
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                            for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                            {
                                                double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                                double zW = zw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))*tan(wings[sD]->alpha(0))/2;
                                                double dxdx2 =-pow(1 + x2_gl_w(jj),-2);
                                                // dydx2 = 0;
                                                double dzdx2 =-tan(wings[sD]->alpha(0))/pow(1 + x2_gl_w(jj), 2);
                                                arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                                double e_11 = dot(e_1, e_1);
                                                double e_12 = dot(e_1, e_2);
                                                double e_22 = dot(e_2, e_2);
                                                double e = e_11*e_22 - pow(e_12, 2);
                                                double e11 = e_22/e;
                                                double e12 =-e_12/e;
                                                double e22 = e_11/e;
                                                for (size_t q = 0; q < wings[sD]->ny; q++)
                                                {
                                                    double t2 = pow(-1, q);
                                                    for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                        for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                        {
                                                            size_t k = i + j*wings[tD]->nx;
                                                            arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                            double r3 = pow(norm(r), 3);
                                                            double dmudxi = dt1 * t2;
                                                            double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                            arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                           dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                          -dydx1(ii)*dxdx2})/sqrt_A;
                                                            arma::vec::fixed<3> J_red = {-n_W(1)*dzdx2,
                                                                                          n_W(0)*dzdx2 - dxdx2*n_W(2),
                                                                                          dxdx2*n_W(1)};
                                                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                            arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                            arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) += dot(q_mu, wings[tD]->nC.row(k)); 
                                                        }
                                                }
                                            }
                                        }
                                    }
                                }
                                else if (c == 1) // East
                                {
                                    arma::vec xw    = T2_gauss * wings[sD]->x.row(wings[sD]->nx-1).t();
                                    arma::vec yw    =-T2_gauss * wings[sD]->y.row(wings[sD]->nx-1).t();
                                    arma::vec zw    = T2_gauss * wings[sD]->z.row(wings[sD]->nx-1).t();
                                    arma::vec dxdx2 = T2_gauss *  dxdx2_gauss.row(wings[sD]->nx-1).t();
                                    arma::vec dydx2 =-T2_gauss *  dydx2_gauss.row(wings[sD]->nx-1).t();
                                    arma::vec dzdx2 = T2_gauss *  dzdx2_gauss.row(wings[sD]->nx-1).t();
                                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                    {
                                        arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                        {
                                            double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                            for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                            {
                                                double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                                double zW = zw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))*tan(wings[sD]->alpha(0))/2;
                                                double dxdx1 = 1/pow(1 - x1_gl_w(ii), 2);
                                                // dydx1 = 0
                                                double dzdx1 = tan(wings[sD]->alpha(0))/pow(1 - x1_gl_w(ii), 2);
                                                arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                                double e_11 = dot(e_1, e_1);
                                                double e_12 = dot(e_1, e_2);
                                                double e_22 = dot(e_2, e_2);
                                                double e = e_11*e_22 - pow(e_12, 2);
                                                double e11 = e_22/e;
                                                double e12 =-e_12/e;
                                                double e22 = e_11/e;
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = dt2;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                        dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                        dxdx1*dydx2(jj)})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                     dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                    -dxdx1*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) += dot(q_mu, wings[tD]->nC.row(k));
                                                    }
                                            }
                                        }
                                    }
                                }
                                else if (c == 2) // North
                                {
                                    arma::vec xw    = T1_gauss * wings[sD]->x.col(wings[sD]->ny-1);
                                    arma::vec yw    =-T1_gauss * wings[sD]->y.col(wings[sD]->ny-1);
                                    arma::vec zw    = T1_gauss * wings[sD]->z.col(wings[sD]->ny-1);
                                    arma::vec dxdx1 = T1_gauss *  dxdx1_gauss.col(wings[sD]->ny-1);
                                    arma::vec dydx1 =-T1_gauss *  dydx1_gauss.col(wings[sD]->ny-1);
                                    arma::vec dzdx1 = T1_gauss *  dzdx1_gauss.col(wings[sD]->ny-1);
                                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                    {
                                        arma::vec::fixed<3> e_1 = {dxdx1(ii), dydx1(ii), dzdx1(ii)};
                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                        {
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                            for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                            {
                                                double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                                double zW = zw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))*tan(wings[sD]->alpha(0))/2;
                                                double dxdx2 = 1/pow(1 - x2_gl_w(jj), 2);
                                                // dydx2 = 0
                                                double dzdx2 = tan(wings[sD]->alpha(0))/pow(1 - x2_gl_w(jj), 2);
                                                arma::vec::fixed<3> e_2 = {dxdx2, 0, dzdx2};
                                                double e_11 = dot(e_1, e_1);
                                                double e_12 = dot(e_1, e_2);
                                                double e_22 = dot(e_2, e_2);
                                                double e = e_11*e_22 - pow(e_12, 2);
                                                double e11 = e_22/e;
                                                double e12 =-e_12/e;
                                                double e22 = e_11/e;
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<3> r = {xW - xC(i, j), yw(ii) - yC(i, j), zW - zC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double dmudxi = dt1;
                                                        double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1(ii), 2) + 2*e12*dzdx1(ii)*dzdx2 + e22*pow(dzdx2, 2)));
                                                        arma::vec::fixed<3> n_W = arma::vec::fixed<3>({dydx1(ii)*dzdx2,
                                                                                                       dzdx1(ii)*dxdx2-dxdx1(ii)*dzdx2,
                                                                                                      -dydx1(ii)*dxdx2})/sqrt_A;
                                                        arma::vec::fixed<3> J_red = {-dzdx2*n_W(1),
                                                                                      dzdx2*n_W(0) - dxdx2*n_W(2),
                                                                                      dxdx2*n_W(1)};
                                                        arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                        arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                        arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) += dot(q_mu, wings[tD]->nC.row(k));
                                                    }
                                            }
                                        }
                                    }
                                }
                                else // West
                                {
                                    arma::vec xw    = T2_gauss * wings[sD]->x.row(0).t();
                                    arma::vec yw    =-T2_gauss * wings[sD]->y.row(0).t();
                                    arma::vec zw    = T2_gauss * wings[sD]->z.row(0).t();
                                    arma::vec dxdx2 = T2_gauss *  dxdx2_gauss.row(0).t();
                                    arma::vec dydx2 =-T2_gauss *  dydx2_gauss.row(0).t();
                                    arma::vec dzdx2 = T2_gauss *  dzdx2_gauss.row(0).t();
                                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                    {
                                        arma::vec::fixed<3> e_2 = {dxdx2(jj), dydx2(jj), dzdx2(jj)};
                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                        {
                                            double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                            for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                            {
                                                double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                                double zW = zw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))*tan(wings[sD]->alpha(0))/2;
                                                double dxdx1 =-1/pow(1 + x1_gl_w(ii), 2);
                                                // dydx1 = 0
                                                double dzdx1 =-tan(wings[sD]->alpha(0))/pow(1 + x1_gl_w(ii), 2);
                                                arma::vec::fixed<3> e_1 = {dxdx1, 0, dzdx1};
                                                double e_11 = dot(e_1, e_1);
                                                double e_12 = dot(e_1, e_2);
                                                double e_22 = dot(e_2, e_2);
                                                double e = e_11*e_22 - pow(e_12, 2);
                                                double e11 = e_22/e;
                                                double e12 =-e_12/e;
                                                double e22 = e_11/e;
                                                for (size_t p = 0; p < wings[sD]->nx; p++)
                                                {
                                                    double t1 = pow(-1, p);
                                                    for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                        for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                        {
                                                            size_t k = i + j*wings[tD]->nx;
                                                            arma::vec::fixed<3> r = {xW - xC(i, j), yw(jj) - yC(i, j), zW - zC(i, j)};
                                                            double r3 = pow(norm(r), 3);
                                                            double dmudxi = t1*dt2;
                                                            double sqrt_A = sqrt(e*(1 + e11*pow(dzdx1, 2) + 2*e12*dzdx1*dzdx2(jj) + e22*pow(dzdx2(jj), 2)));
                                                            arma::vec::fixed<3> n_W = arma::vec::fixed<3>({-dzdx1*dydx2(jj),
                                                                                                            dzdx1*dxdx2(jj)-dxdx1*dzdx2(jj),
                                                                                                            dxdx1*dydx2(jj)})/sqrt_A;
                                                            arma::vec::fixed<3> J_red = {n_W(1)*dzdx1,
                                                                                         dxdx1*n_W(2) - n_W(0)*dzdx1,
                                                                                        -dxdx1*n_W(1)};
                                                            arma::vec::fixed<3> gradmu = J_red*dmudxi;
                                                            arma::vec::fixed<3> gamma_W = cross(gradmu, n_W);
                                                            arma::vec q_mu = gl_x_w[ii].weight * gl_y_w[jj].weight * cross(gamma_W, r)/r3;
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx, 0) += dot(q_mu, wings[tD]->nC.row(k)); 
                                                        }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
    }
    solve();
}