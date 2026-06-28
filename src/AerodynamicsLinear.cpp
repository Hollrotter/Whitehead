#include "Aerodynamics.hpp"

void Aerodynamics::linear()
{
    analysis = Analysis::linear;
    // Influence of the wing surfaces on each other
    bw.set_size(wings.size(), wings.size());
    for (size_t sD = 0; sD < wings.size(); sD++)
    {
        wings[sD]->analysis = Analysis::linear;
        for (size_t tD = 0; tD < wings.size(); tD++)
            if (sD != tD)
            {
                std::vector<fastgl::QuadPair> gl_x(wings[sD]->nx), gl_y(wings[sD]->ny);
                arma::vec x1_gl(wings[sD]->nx, arma::fill::none), x2_gl(wings[sD]->ny, arma::fill::none);
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
                auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                arma::mat xC = wings[tD]->xC;
                arma::mat yC = wings[tD]->yC;
                bw(tD, sD).set_size(wings[tD]->nxy, wings[sD]->nxy);
                for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                    for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                    {
                        size_t k = i+j*wings[tD]->nx;
                        for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                            for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                            {
                                arma::vec::fixed<2> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j)};
                                double r3 = pow(norm(r), 3);
                                double t2   = 1; // First Chebyshev Polynomial (j)
                                double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

                                double dt2   = 0;
                                double dt2p1 = 1;
                                for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction of source
                                {
                                    double t1   = 1; // First Chebyshev Polynomial (i)
                                    double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                                    double dt1   = 0;
                                    double dt1p1 = 1;
                                    for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction of source
                                    {
                                        double dmudx1 = dt1 *  t2;
                                        double dmudx2 =  t1 * dt2;
                                        bw(tD, sD)(k, p+q*wings[sD]->nx) -= gl_x[ii].weight * gl_y[jj].weight / r3
                                            *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                            - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                                        std::swap(t1, t1p1); // Swap order so the oldes Chebyshev Polynomial will be overwritten
                                        t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1); // Calculate next Chebyshev Polynomial

                                        std::swap(dt1, dt1p1);
                                        dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                                    }
                                    std::swap(t2, t2p1); // Swap order so the oldest Chebyshev Polynomial will be overwritten
                                    t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1); // Calculate next Chebyshev Polynomial

                                    std::swap(dt2, dt2p1);
                                    dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
                                }
                            }
                    }
                for (const Wake* w:wings[sD]->wakes)
                    for (size_t c = 0; c < 4; c++)
                        if (wings[sD]->chi[c] == w->chi)
                        {
                            std::vector<fastgl::QuadPair> gl_x_w(wings[sD]->nx), gl_y_w(wings[sD]->ny);
                            arma::vec x1_gl_w(wings[sD]->nx, arma::fill::none), x2_gl_w(wings[sD]->ny, arma::fill::none);
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
                            arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                            arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                            if (c == 0) // South
                            {
                                arma::vec xw = T1_gauss * wings[sD]->x.col(0);
                                arma::vec yw = T1_gauss * wings[sD]->y.col(0);
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                            for (size_t q = 0; q < wings[sD]->ny; q++)
                                            {
                                                double t2 = pow(-1, q);
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx) -= gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*t2*r(1)/pow(1 + x2_gl_w(jj), 2)/r3; 
                                                    }
                                            }
                                        }
                                    }
                            }
                            else if (c == 1) // East
                            {
                                arma::vec xw = T2_gauss * wings[sD]->x.row(wings[sD]->nx-1).t();
                                arma::vec yw = T2_gauss * wings[sD]->y.row(wings[sD]->nx-1).t();
                                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                            for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                {
                                                    size_t k = i + j*wings[tD]->nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double I = gl_x_w[ii].weight*gl_y_w[jj].weight*dt2*r(1)/pow(1 - x1_gl_w(ii), 2)/r3;
                                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx) -= I;
                                                }
                                        }
                                    }
                            }
                            else if (c == 2) // North
                            {
                                arma::vec xw = T1_gauss * wings[sD]->x.col(wings[sD]->ny-1);
                                arma::vec yw = T1_gauss * wings[sD]->y.col(wings[sD]->ny-1);
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                    for (size_t p = 0; p < wings[sD]->nx; p++)
                                    {
                                        double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                        for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        {
                                            double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                            for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                {
                                                    size_t k = i + j*wings[tD]->nx;
                                                    arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                    double r3 = pow(norm(r), 3);
                                                    double I =-gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*r(1)/pow(1 - x2_gl_w(jj), 2)/r3;
                                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx) -= I;
                                                }
                                        }
                                    }
                            }
                            else // West
                            {
                                arma::vec xw = T2_gauss * wings[sD]->x.row(0).t();
                                arma::vec yw = T2_gauss * wings[sD]->y.row(0).t();
                                for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                    for (size_t q = 0; q < wings[sD]->ny; q++)
                                    {
                                        double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                        for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        {
                                            double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                            for (size_t p = 0; p < wings[sD]->nx; p++)
                                            {
                                                double t1 = pow(-1, p);
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        bw(tD, sD)(k, p+q*wings[sD]->nx) += gl_x_w[ii].weight*gl_y_w[jj].weight*t1*dt2*r(1)/pow(1 + x1_gl_w(ii), 2)/r3;
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
                    arma::vec x1_gl(wings[sD]->nx, arma::fill::none), x2_gl(wings[sD]->ny, arma::fill::none);
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
                    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, wings[sD]->chi);
                    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, wings[sD]->chi);
                    y_gl     =-y_gl;
                    dy_gldx1 =-dy_gldx1;
                    dy_gldx2 =-dy_gldx2;
                    arma::mat xC = wings[tD]->xC;
                    arma::mat yC = wings[tD]->yC;
                    for (size_t j = 1; j < wings[tD]->ny-1; j++) // Loop over Collocation Points in 2-direction of target
                        for (size_t i = 1; i < wings[tD]->nx-1; i++) // Loop over Collocation Points in 1-direction of target
                        {
                            size_t k = i+j*wings[tD]->nx;
                            for (size_t jj = 0; jj < wings[sD]->ny; jj++) // Loop over Legendre nodes 2-direction of source
                                for (size_t ii = 0; ii < wings[sD]->nx; ii++) // Loop over Legendre nodes 1-direction of source
                                {
                                    arma::vec::fixed<2> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j)};
                                    double r3 = pow(norm(r), 3);
                                    double t2   = 1; // First Chebyshev Polynomial (j)
                                    double t2p1 = x2_gl(jj); // Second Chebyshev Polynomial (j+1 -> jp1)

                                    double dt2   = 0;
                                    double dt2p1 = 1;
                                    for (size_t q = 0; q < wings[sD]->ny; q++) // Loop over Chebyshev Polynomial 2-direction of source
                                    {
                                        double t1   = 1; // First Chebyshev Polynomial (i)
                                        double t1p1 = x1_gl(ii); // Second Chebyshev Polynomial (i+1 -> ip1)

                                        double dt1   = 0;
                                        double dt1p1 = 1;
                                        for (size_t p = 0; p < wings[sD]->nx; p++) // Loop over Chebyshev Polynomial 1-direction of source
                                        {
                                            double dmudx1 = dt1 *  t2;
                                            double dmudx2 =  t1 * dt2;
                                            bw(tD, sD)(k, p+q*wings[sD]->nx) += gl_x[ii].weight * gl_y[jj].weight / r3
                                                *(r(0)*(dy_gldx2(ii, jj)*dmudx1 - dy_gldx1(ii, jj)*dmudx2)
                                                - r(1)*(dx_gldx2(ii, jj)*dmudx1 - dx_gldx1(ii, jj)*dmudx2));
                                            std::swap(t1, t1p1); // Swap order so the oldes Chebyshev Polynomial will be overwritten
                                            t1p1 = boost::math::chebyshev_next(x1_gl(ii), t1, t1p1); // Calculate next Chebyshev Polynomial

                                            std::swap(dt1, dt1p1);
                                            dt1p1 = (p == 0) ? 4*x1_gl(ii) : (p+2)*(2*t1 + dt1p1/p);
                                        }
                                        std::swap(t2, t2p1); // Swap order so the oldest Chebyshev Polynomial will be overwritten
                                        t2p1 = boost::math::chebyshev_next(x2_gl(jj), t2, t2p1); // Calculate next Chebyshev Polynomial

                                        std::swap(dt2, dt2p1);
                                        dt2p1 = (q == 0) ? 4*x2_gl(jj) : (q+2)*(2*t2 + dt2p1/q);
                                    }
                                }
                        }
                    for (const Wake* w:wings[sD]->wakes)
                        for (size_t c = 0; c < 4; c++)
                            if (wings[sD]->chi[c] == w->chi)
                            {
                                std::vector<fastgl::QuadPair> gl_x_w(wings[sD]->nx), gl_y_w(wings[sD]->ny);
                                arma::vec x1_gl_w(wings[sD]->nx, arma::fill::none), x2_gl_w(wings[sD]->ny, arma::fill::none);
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
                                arma::mat T1_gauss = Lagrange::interpolationMatrix(x1_gauss, x1_gl_w);
                                arma::mat T2_gauss = Lagrange::interpolationMatrix(x2_gauss, x2_gl_w);
                                if (c == 0) // South
                                {
                                    arma::vec xw = T1_gauss * wings[sD]->x.col(0);
                                    arma::vec yw =-T1_gauss * wings[sD]->y.col(0);
                                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                        {
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                            for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                            {
                                                double xW = xw(ii) + (1 - x2_gl_w(jj))/(1 + x2_gl_w(jj))/2;
                                                for (size_t q = 0; q < wings[sD]->ny; q++)
                                                {
                                                    double t2 = pow(-1, q);
                                                    for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                        for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                        {
                                                            size_t k = i + j*wings[tD]->nx;
                                                            arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                            double r3 = pow(norm(r), 3);
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx) += gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*t2*r(1)/pow(1 + x2_gl_w(jj), 2)/r3; 
                                                        }
                                                }
                                            }
                                        }
                                }
                                else if (c == 1) // East
                                {
                                    arma::vec xw = T2_gauss * wings[sD]->x.row(wings[sD]->nx-1).t();
                                    arma::vec yw =-T2_gauss * wings[sD]->y.row(wings[sD]->nx-1).t();
                                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                        {
                                            double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                            for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                            {
                                                double xW = xw(jj) + (1 + x1_gl_w(ii))/(1 - x1_gl_w(ii))/2;
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double I = gl_x_w[ii].weight*gl_y_w[jj].weight*dt2*r(1)/pow(1 - x1_gl_w(ii), 2)/r3;
                                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx) += I;
                                                    }
                                            }
                                        }
                                }
                                else if (c == 2) // North
                                {
                                    arma::vec xw = T1_gauss * wings[sD]->x.col(wings[sD]->ny-1);
                                    arma::vec yw =-T1_gauss * wings[sD]->y.col(wings[sD]->ny-1);
                                    for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                        for (size_t p = 0; p < wings[sD]->nx; p++)
                                        {
                                            double dt1 = boost::math::chebyshev_t_prime(p, x1_gl_w(ii));
                                            for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                            {
                                                double xW = xw(ii) + (1 + x2_gl_w(jj))/(1 - x2_gl_w(jj))/2;
                                                for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                    for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                    {
                                                        size_t k = i + j*wings[tD]->nx;
                                                        arma::vec::fixed<2> r = {xW - xC(i, j), yw(ii) - yC(i, j)};
                                                        double r3 = pow(norm(r), 3);
                                                        double I =-gl_x_w[ii].weight*gl_y_w[jj].weight*dt1*r(1)/pow(1 - x2_gl_w(jj), 2)/r3;
                                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx) += I;
                                                    }
                                            }
                                        }
                                }
                                else // West
                                {
                                arma::vec xw = T2_gauss * wings[sD]->x.row(0).t();
                                arma::vec yw =-T2_gauss * wings[sD]->y.row(0).t();
                                    for (size_t jj = 0; jj < wings[sD]->ny; jj++)
                                        for (size_t q = 0; q < wings[sD]->ny; q++)
                                        {
                                            double dt2 = boost::math::chebyshev_t_prime(q, x2_gl_w(jj));
                                            for (size_t ii = 0; ii < wings[sD]->nx; ii++)
                                            {
                                                double xW = xw(jj) + (1 - x1_gl_w(ii))/(1 + x1_gl_w(ii))/2;
                                                for (size_t p = 0; p < wings[sD]->nx; p++)
                                                {
                                                    double t1 = pow(-1, p);
                                                    for (size_t j = 1; j < wings[tD]->ny-1; j++)
                                                        for (size_t i = 1; i < wings[tD]->nx-1; i++)
                                                        {
                                                            size_t k = i + j*wings[tD]->nx;
                                                            arma::vec::fixed<2> r = {xW - xC(i, j), yw(jj) - yC(i, j)};
                                                            double r3 = pow(norm(r), 3);
                                                            bw(tD, sD)(k, p+q*wings[sD]->nx) -= gl_x_w[ii].weight*gl_y_w[jj].weight*t1*dt2*r(1)/pow(1 + x1_gl_w(ii), 2)/r3;
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