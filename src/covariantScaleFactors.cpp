#include "Lagrange.hpp"

std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
    Lagrange::covariantScaleFactors(const std::array<CurveInterpolant*, 4> chi)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactors(x1, x2, chi);
}

std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
    Lagrange::covariantScaleFactors(const arma::vec &x1, const arma::vec &x2, const std::array<CurveInterpolant*, 4> chi)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactors", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::rowvec dxdx1_east(x2.size(), arma::fill::none), dydx1_east(x2.size(), arma::fill::none),
                 dxdx2_east(x2.size(), arma::fill::none), dydx2_east(x2.size(), arma::fill::none),
                 dxdx1_west(x2.size(), arma::fill::none), dydx1_west(x2.size(), arma::fill::none),
                 dxdx2_west(x2.size(), arma::fill::none), dydx2_west(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto [X_1,  Y_1]  = chi[0]->evaluate(-sign1);
        auto [X_3,  Y_3]  = chi[2]->evaluate(-sign3);
        auto [Xs_1, Ys_1] = chi[0]->derivative(-sign1);
        auto [Xs_3, Ys_3] = chi[2]->derivative(-sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_4, Ys_4] = chi[3]->derivative(sign4*x2j);
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        dxdx1_west(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_west(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_west(j) = (2*Xs_4 + X_3 - X_1)/2 - (x_4 - x_1)/2;
        dydx2_west(j) = (2*Ys_4 + Y_3 - Y_1)/2 - (y_4 - y_1)/2;

        std::tie(X_1,  Y_1)  = chi[0]->evaluate(sign1);
        std::tie(X_3,  Y_3)  = chi[2]->evaluate(sign3);
        std::tie(Xs_1, Ys_1) = chi[0]->derivative(sign1);
        std::tie(Xs_3, Ys_3) = chi[2]->derivative(sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        x2j = x2(j);
        std::tie(X_2,  Y_2)  = chi[1]->evaluate(sign2*x2j);
        std::tie(X_4,  Y_4)  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_2, Ys_2] = chi[1]->derivative(sign2*x2j);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        dxdx1_east(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_east(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_east(j) = (2*Xs_2 + X_3 - X_1)/2 - (x_3 - x_2)/2;
        dydx2_east(j) = (2*Ys_2 + Y_3 - Y_1)/2 - (y_3 - y_2)/2;
    }
    arma::vec dxdx1_south(x1.size(), arma::fill::none), dydx1_south(x1.size(), arma::fill::none),
              dxdx2_south(x1.size(), arma::fill::none), dydx2_south(x1.size(), arma::fill::none),
              dxdx1_north(x1.size(), arma::fill::none), dydx1_north(x1.size(), arma::fill::none),
              dxdx2_north(x1.size(), arma::fill::none), dydx2_north(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto [X_2,  Y_2]  = chi[1]->evaluate(-sign2);
        auto [X_4,  Y_4]  = chi[3]->evaluate(-sign4);
        auto [Xs_2, Ys_2] = chi[1]->derivative(-sign2);
        auto [Xs_4, Ys_4] = chi[3]->derivative(-sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_1, Ys_1] = chi[0]->derivative(sign1*x1i);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        dxdx1_south(i) = (X_2 - X_4 + 2*Xs_1)/2 - (x_2 - x_1)/2;
        dydx1_south(i) = (Y_2 - Y_4 + 2*Ys_1)/2 - (y_2 - y_1)/2;
        dxdx2_south(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_south(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;

        std::tie(X_2,  Y_2)  = chi[1]->evaluate(sign2);
        std::tie(X_4,  Y_4)  = chi[3]->evaluate(sign4);
        std::tie(Xs_2, Ys_2) = chi[1]->derivative(sign2);
        std::tie(Xs_4, Ys_4) = chi[3]->derivative(sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        x1i = x1(i);
        std::tie(X_1,  Y_1)  = chi[0]->evaluate(sign1*x1i);
        std::tie(X_3,  Y_3)  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_3, Ys_3] = chi[2]->derivative(sign3*x1i);
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        dxdx1_north(i) = (X_2 - X_4 + 2*Xs_3)/2 - (x_3 - x_4)/2;
        dydx1_north(i) = (Y_2 - Y_4 + 2*Ys_3)/2 - (y_3 - y_4)/2;
        dxdx2_north(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_north(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }
    arma::vec e_11_south = pow(dxdx1_south, 2) + pow(dydx1_south, 2);
    arma::vec e_12_south = dxdx1_south%dxdx2_south + dydx1_south%dydx2_south;
    arma::vec e_22_south = pow(dxdx2_south, 2) + pow(dydx2_south, 2);

    arma::vec e_11_north = pow(dxdx1_north, 2) + pow(dydx1_north, 2);
    arma::vec e_12_north = dxdx1_north%dxdx2_north + dydx1_north%dydx2_north;
    arma::vec e_22_north = pow(dxdx2_north, 2) + pow(dydx2_north, 2);

    arma::rowvec e_11_west = pow(dxdx1_west, 2) + pow(dydx1_west, 2);
    arma::rowvec e_12_west = dxdx1_west%dxdx2_west + dydx1_west%dydx2_west;
    arma::rowvec e_22_west = pow(dxdx2_west, 2) + pow(dydx2_west, 2);

    arma::rowvec e_11_east = pow(dxdx1_east, 2) + pow(dydx1_east, 2);
    arma::rowvec e_12_east = dxdx1_east%dxdx2_east + dydx1_east%dydx2_east;
    arma::rowvec e_22_east = pow(dxdx2_east, 2) + pow(dydx2_east, 2);

    arma::vec    e_south = e_11_south%e_22_south - pow(e_12_south, 2);
    arma::vec    e_north = e_11_north%e_22_north - pow(e_12_north, 2);
    arma::rowvec e_west  = e_11_west %e_22_west  - pow(e_12_west,  2);
    arma::rowvec e_east  = e_11_east %e_22_east  - pow(e_12_east,  2);

    arma::vec e12_south =-e_12_south/e_south;
    arma::vec e22_south = e_11_south/e_south;

    arma::vec e12_north =-e_12_north/e_north;
    arma::vec e22_north = e_11_north/e_north;

    arma::rowvec e11_west = e_22_west/e_west;
    arma::rowvec e12_west =-e_12_west/e_west;

    arma::rowvec e11_east = e_22_east/e_east;
    arma::rowvec e12_east =-e_12_east/e_east;

    arma::rowvec h_1s1_east = sqrt(e11_east);
    arma::rowvec h_1s2_east = e12_east/sqrt(e11_east);
    arma::rowvec h_1s1_west = sqrt(e11_west);
    arma::rowvec h_1s2_west = e12_west/sqrt(e11_west);
    arma::vec h_2s2_south = sqrt(e22_south);
    arma::vec h_2s1_south = e12_south/sqrt(e22_south);
    arma::vec h_2s2_north = sqrt(e22_north);
    arma::vec h_2s1_north = e12_north/sqrt(e22_north);

    return std::make_tuple(h_2s2_south, h_2s1_south, h_1s1_east, h_1s2_east, h_2s2_north, h_2s1_north, h_1s1_west, h_1s2_west);
}

std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
    Lagrange::covariantScaleFactors(const std::array<CurveInterpolant*, 4> chi, const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactors(x1, x2, chi, dzdx1, dzdx2);
}

std::tuple<arma::vec, arma::vec, arma::rowvec, arma::rowvec, arma::vec, arma::vec, arma::rowvec, arma::rowvec>
    Lagrange::covariantScaleFactors(const arma::vec &x1, const arma::vec &x2, const std::array<CurveInterpolant*, 4> chi,
        const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactors", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);
    arma::rowvec dxdx1_east(x2.size(), arma::fill::none), dydx1_east(x2.size(), arma::fill::none),
                 dxdx2_east(x2.size(), arma::fill::none), dydx2_east(x2.size(), arma::fill::none),
                 dxdx1_west(x2.size(), arma::fill::none), dydx1_west(x2.size(), arma::fill::none),
                 dxdx2_west(x2.size(), arma::fill::none), dydx2_west(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto [X_1,  Y_1]  = chi[0]->evaluate(-sign1);
        auto [X_3,  Y_3]  = chi[2]->evaluate(-sign3);
        auto [Xs_1, Ys_1] = chi[0]->derivative(-sign1);
        auto [Xs_3, Ys_3] = chi[2]->derivative(-sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_4, Ys_4] = chi[3]->derivative(sign4*x2j);
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        dxdx1_west(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_west(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_west(j) = (2*Xs_4 + X_3 - X_1)/2 - (x_4 - x_1)/2;
        dydx2_west(j) = (2*Ys_4 + Y_3 - Y_1)/2 - (y_4 - y_1)/2;

        std::tie(X_1,  Y_1)  = chi[0]->evaluate(sign1);
        std::tie(X_3,  Y_3)  = chi[2]->evaluate(sign3);
        std::tie(Xs_1, Ys_1) = chi[0]->derivative(sign1);
        std::tie(Xs_3, Ys_3) = chi[2]->derivative(sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        x2j = x2(j);
        std::tie(X_2,  Y_2)  = chi[1]->evaluate(sign2*x2j);
        std::tie(X_4,  Y_4)  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_2, Ys_2] = chi[1]->derivative(sign2*x2j);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        dxdx1_east(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_east(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_east(j) = (2*Xs_2 + X_3 - X_1)/2 - (x_3 - x_2)/2;
        dydx2_east(j) = (2*Ys_2 + Y_3 - Y_1)/2 - (y_3 - y_2)/2;
    }
    arma::vec dxdx1_south(x1.size(), arma::fill::none), dydx1_south(x1.size(), arma::fill::none),
              dxdx2_south(x1.size(), arma::fill::none), dydx2_south(x1.size(), arma::fill::none),
              dxdx1_north(x1.size(), arma::fill::none), dydx1_north(x1.size(), arma::fill::none),
              dxdx2_north(x1.size(), arma::fill::none), dydx2_north(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto [X_2,  Y_2]  = chi[1]->evaluate(-sign2);
        auto [X_4,  Y_4]  = chi[3]->evaluate(-sign4);
        auto [Xs_2, Ys_2] = chi[1]->derivative(-sign2);
        auto [Xs_4, Ys_4] = chi[3]->derivative(-sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_1, Ys_1] = chi[0]->derivative(sign1*x1i);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        dxdx1_south(i) = (X_2 - X_4 + 2*Xs_1)/2 - (x_2 - x_1)/2;
        dydx1_south(i) = (Y_2 - Y_4 + 2*Ys_1)/2 - (y_2 - y_1)/2;
        dxdx2_south(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_south(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;

        std::tie(X_2,  Y_2)  = chi[1]->evaluate(sign2);
        std::tie(X_4,  Y_4)  = chi[3]->evaluate(sign4);
        std::tie(Xs_2, Ys_2) = chi[1]->derivative(sign2);
        std::tie(Xs_4, Ys_4) = chi[3]->derivative(sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        x1i = x1(i);
        std::tie(X_1,  Y_1)  = chi[0]->evaluate(sign1*x1i);
        std::tie(X_3,  Y_3)  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_3, Ys_3] = chi[2]->derivative(sign3*x1i);
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        dxdx1_north(i) = (X_2 - X_4 + 2*Xs_3)/2 - (x_3 - x_4)/2;
        dydx1_north(i) = (Y_2 - Y_4 + 2*Ys_3)/2 - (y_3 - y_4)/2;
        dxdx2_north(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_north(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }

    arma::vec    a_11_south = pow(dxdx1_south, 2) + pow(dydx1_south, 2) + pow(dzdx1.col(0), 2);
    arma::vec    a_12_south = dxdx1_south%dxdx2_south + dydx1_south%dydx2_south + dzdx1.col(0)%dzdx2.col(0);
    arma::vec    a_22_south = pow(dxdx2_south, 2) + pow(dydx2_south, 2) + pow(dzdx2.col(0), 2);

    arma::vec    a_11_north = pow(dxdx1_north, 2) + pow(dydx1_north, 2) + pow(dzdx1.col(x2.size()-1), 2);
    arma::vec    a_12_north = dxdx1_north%dxdx2_north + dydx1_north%dydx2_north + dzdx1.col(x2.size()-1)%dzdx2.col(x2.size()-1);
    arma::vec    a_22_north = pow(dxdx2_north, 2) + pow(dydx2_north, 2) + pow(dzdx2.col(x2.size()-1), 2);

    arma::rowvec a_11_west = pow(dxdx1_west, 2) + pow(dydx1_west, 2) + pow(dzdx1.row(0), 2);
    arma::rowvec a_12_west = dxdx1_west%dxdx2_west + dydx1_west%dydx2_west + dzdx1.row(0)%dzdx2.row(0);
    arma::rowvec a_22_west = pow(dxdx2_west, 2) + pow(dydx2_west, 2) + pow(dzdx2.row(0), 2);

    arma::rowvec a_11_east = pow(dxdx1_east, 2) + pow(dydx1_east, 2) + pow(dzdx1.row(x1.size()-1), 2);
    arma::rowvec a_12_east = dxdx1_east%dxdx2_east + dydx1_east%dydx2_east + dzdx1.row(x1.size()-1)%dzdx2.row(x1.size()-1);
    arma::rowvec a_22_east = pow(dxdx2_east, 2) + pow(dydx2_east, 2) + pow(dzdx2.row(x1.size()-1), 2);

    arma::vec    a_south = a_11_south%a_22_south - pow(a_12_south, 2);
    arma::vec    a_north = a_11_north%a_22_north - pow(a_12_north, 2);
    arma::rowvec a_west  = a_11_west %a_22_west  - pow(a_12_west,  2);
    arma::rowvec a_east  = a_11_east %a_22_east  - pow(a_12_east,  2);

    arma::vec    a12_south =-a_12_south/a_south;
    arma::vec    a22_south = a_11_south/a_south;

    arma::vec    a12_north =-a_12_north/a_north;
    arma::vec    a22_north = a_11_north/a_north;

    arma::rowvec a11_west = a_22_west/a_west;
    arma::rowvec a12_west =-a_12_west/a_west;

    arma::rowvec a11_east = a_22_east/a_east;
    arma::rowvec a12_east =-a_12_east/a_east;

    arma::rowvec h_1s1_east  = sqrt(a11_east);
    arma::rowvec h_1s2_east  = a12_east/sqrt(a11_east);
    arma::rowvec h_1s1_west  = sqrt(a11_west);
    arma::rowvec h_1s2_west  = a12_west/sqrt(a11_west);
    arma::vec    h_2s2_south = sqrt(a22_south);
    arma::vec    h_2s1_south = a12_south/sqrt(a22_south);
    arma::vec    h_2s2_north = sqrt(a22_north);
    arma::vec    h_2s1_north = a12_north/sqrt(a22_north);

    return std::make_tuple(h_2s2_south, h_2s1_south, h_1s1_east, h_1s2_east, h_2s2_north, h_2s1_north, h_1s1_west, h_1s2_west);
}

std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsSouth(const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    return covariantScaleFactorsSouth(x1, chi);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsSouth(const arma::vec &x1, const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsSouth", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::vec dxdx1_south(x1.size(), arma::fill::none), dydx1_south(x1.size(), arma::fill::none),
              dxdx2_south(x1.size(), arma::fill::none), dydx2_south(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto [X_2,  Y_2]  = chi[1]->evaluate(-sign2);
        auto [X_4,  Y_4]  = chi[3]->evaluate(-sign4);
        auto [Xs_2, Ys_2] = chi[1]->derivative(-sign2);
        auto [Xs_4, Ys_4] = chi[3]->derivative(-sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_1, Ys_1] = chi[0]->derivative(sign1*x1i);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        dxdx1_south(i) = (X_2 - X_4 + 2*Xs_1)/2 - (x_2 - x_1)/2;
        dydx1_south(i) = (Y_2 - Y_4 + 2*Ys_1)/2 - (y_2 - y_1)/2;
        dxdx2_south(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_south(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }
    arma::vec e_11_south = pow(dxdx1_south, 2) + pow(dydx1_south, 2);
    arma::vec e_12_south = dxdx1_south%dxdx2_south + dydx1_south%dydx2_south;
    arma::vec e_22_south = pow(dxdx2_south, 2) + pow(dydx2_south, 2);

    arma::vec e_south = e_11_south%e_22_south - pow(e_12_south, 2);

    arma::vec e12_south =-e_12_south/e_south;
    arma::vec e22_south = e_11_south/e_south;

    arma::vec h_2s2_south = sqrt(e22_south);
    arma::vec h_2s1_south = e12_south/sqrt(e22_south);

    return std::make_tuple(h_2s2_south, h_2s1_south);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsSouth(const std::array<Lagrange::CurveInterpolant*, 4> chi, const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    return covariantScaleFactorsSouth(x1, chi, dzdx1, dzdx2);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsSouth(const arma::vec &x1, const std::array<Lagrange::CurveInterpolant*, 4> chi,
    const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsSouth", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::vec dxdx1_south(x1.size(), arma::fill::none), dydx1_south(x1.size(), arma::fill::none),
              dxdx2_south(x1.size(), arma::fill::none), dydx2_south(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto [X_2,  Y_2]  = chi[1]->evaluate(-sign2);
        auto [X_4,  Y_4]  = chi[3]->evaluate(-sign4);
        auto [Xs_2, Ys_2] = chi[1]->derivative(-sign2);
        auto [Xs_4, Ys_4] = chi[3]->derivative(-sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_1, Ys_1] = chi[0]->derivative(sign1*x1i);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        dxdx1_south(i) = (X_2 - X_4 + 2*Xs_1)/2 - (x_2 - x_1)/2;
        dydx1_south(i) = (Y_2 - Y_4 + 2*Ys_1)/2 - (y_2 - y_1)/2;
        dxdx2_south(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_south(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }

    arma::vec a_11_south = pow(dxdx1_south, 2) + pow(dydx1_south, 2) + pow(dzdx1.col(0), 2);
    arma::vec a_12_south = dxdx1_south%dxdx2_south + dydx1_south%dydx2_south + dzdx1.col(0)%dzdx2.col(0);
    arma::vec a_22_south = pow(dxdx2_south, 2) + pow(dydx2_south, 2) + pow(dzdx2.col(0), 2);

    arma::vec a_south = a_11_south%a_22_south - pow(a_12_south, 2);

    arma::vec a12_south =-a_12_south/a_south;
    arma::vec a22_south = a_11_south/a_south;

    arma::vec h_2s2_south = sqrt(a22_south);
    arma::vec h_2s1_south = a12_south/sqrt(a22_south);

    return std::make_tuple(h_2s2_south, h_2s1_south);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsEast(const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactorsEast(x2, chi);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsEast(const arma::vec &x2, const std::array<CurveInterpolant*, 4> chi)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsEast", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::rowvec dxdx1_east(x2.size(), arma::fill::none), dydx1_east(x2.size(), arma::fill::none),
                 dxdx2_east(x2.size(), arma::fill::none), dydx2_east(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto[X_1,  Y_1]  = chi[0]->evaluate(sign1);
        auto[X_3,  Y_3]  = chi[2]->evaluate(sign3);
        auto[Xs_1, Ys_1] = chi[0]->derivative(sign1);
        auto[Xs_3, Ys_3] = chi[2]->derivative(sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_2, Ys_2] = chi[1]->derivative(sign2*x2j);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        dxdx1_east(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_east(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_east(j) = (2*Xs_2 + X_3 - X_1)/2 - (x_3 - x_2)/2;
        dydx2_east(j) = (2*Ys_2 + Y_3 - Y_1)/2 - (y_3 - y_2)/2;
    }

    arma::rowvec e_11_east = pow(dxdx1_east, 2) + pow(dydx1_east, 2);
    arma::rowvec e_12_east = dxdx1_east%dxdx2_east + dydx1_east%dydx2_east;
    arma::rowvec e_22_east = pow(dxdx2_east, 2) + pow(dydx2_east, 2);

    arma::rowvec e_east = e_11_east %e_22_east  - pow(e_12_east,  2);

    arma::rowvec e11_east = e_22_east/e_east;
    arma::rowvec e12_east =-e_12_east/e_east;

    arma::rowvec h_1s1_east = sqrt(e11_east);
    arma::rowvec h_1s2_east = e12_east/sqrt(e11_east);

    return std::make_tuple(h_1s1_east, h_1s2_east);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsEast(const std::array<Lagrange::CurveInterpolant*, 4> chi, const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactorsEast(x2, chi, dzdx1, dzdx2);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsEast(const arma::vec &x2, const std::array<Lagrange::CurveInterpolant*, 4> chi,
    const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsEast", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);
    arma::rowvec dxdx1_east(x2.size(), arma::fill::none), dydx1_east(x2.size(), arma::fill::none),
                 dxdx2_east(x2.size(), arma::fill::none), dydx2_east(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto[X_1,  Y_1]  = chi[0]->evaluate(sign1);
        auto[X_3,  Y_3]  = chi[2]->evaluate(sign3);
        auto[Xs_1, Ys_1] = chi[0]->derivative(sign1);
        auto[Xs_3, Ys_3] = chi[2]->derivative(sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_2, Ys_2] = chi[1]->derivative(sign2*x2j);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        dxdx1_east(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_east(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_east(j) = (2*Xs_2 + X_3 - X_1)/2 - (x_3 - x_2)/2;
        dydx2_east(j) = (2*Ys_2 + Y_3 - Y_1)/2 - (y_3 - y_2)/2;
    }

    arma::rowvec a_11_east = pow(dxdx1_east, 2) + pow(dydx1_east, 2) + pow(dzdx1.row(dzdx1.n_rows-1), 2);
    arma::rowvec a_12_east = dxdx1_east%dxdx2_east + dydx1_east%dydx2_east + dzdx1.row(dzdx1.n_rows-1)%dzdx2.row(dzdx2.n_rows-1);
    arma::rowvec a_22_east = pow(dxdx2_east, 2) + pow(dydx2_east, 2) + pow(dzdx2.row(dzdx2.n_rows-1), 2);

    arma::rowvec a_east = a_11_east %a_22_east - pow(a_12_east,  2);

    arma::rowvec a11_east = a_22_east/a_east;
    arma::rowvec a12_east =-a_12_east/a_east;

    arma::rowvec h_1s1_east  = sqrt(a11_east);
    arma::rowvec h_1s2_east  = a12_east/sqrt(a11_east);

    return std::make_tuple(h_1s1_east, h_1s2_east);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsNorth(const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    return covariantScaleFactorsNorth(x1, chi);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsNorth(const arma::vec &x1, const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsNorth", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::vec dxdx1_north(x1.size(), arma::fill::none), dydx1_north(x1.size(), arma::fill::none),
              dxdx2_north(x1.size(), arma::fill::none), dydx2_north(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto[X_2,  Y_2]  = chi[1]->evaluate(sign2);
        auto[X_4,  Y_4]  = chi[3]->evaluate(sign4);
        auto[Xs_2, Ys_2] = chi[1]->derivative(sign2);
        auto[Xs_4, Ys_4] = chi[3]->derivative(sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_3, Ys_3] = chi[2]->derivative(sign3*x1i);
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        dxdx1_north(i) = (X_2 - X_4 + 2*Xs_3)/2 - (x_3 - x_4)/2;
        dydx1_north(i) = (Y_2 - Y_4 + 2*Ys_3)/2 - (y_3 - y_4)/2;
        dxdx2_north(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_north(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }

    arma::vec e_11_north = pow(dxdx1_north, 2) + pow(dydx1_north, 2);
    arma::vec e_12_north = dxdx1_north%dxdx2_north + dydx1_north%dydx2_north;
    arma::vec e_22_north = pow(dxdx2_north, 2) + pow(dydx2_north, 2);

    arma::vec    e_north = e_11_north%e_22_north - pow(e_12_north, 2);

    arma::vec e12_north =-e_12_north/e_north;
    arma::vec e22_north = e_11_north/e_north;

    arma::vec h_2s2_north = sqrt(e22_north);
    arma::vec h_2s1_north = e12_north/sqrt(e22_north);

    return std::make_tuple(h_2s2_north, h_2s1_north);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsNorth(const std::array<Lagrange::CurveInterpolant*, 4> chi, const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    arma::vec x1 = Chebyshev::gaussLobatto(chi[0]->getNodes().size());
    return covariantScaleFactorsNorth(x1, chi, dzdx1, dzdx2);
}
std::tuple<arma::vec, arma::vec> Lagrange::covariantScaleFactorsNorth(const arma::vec &x1, const std::array<Lagrange::CurveInterpolant*, 4> chi,
    const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsNorth", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::vec dxdx1_north(x1.size(), arma::fill::none), dydx1_north(x1.size(), arma::fill::none),
              dxdx2_north(x1.size(), arma::fill::none), dydx2_north(x1.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < x1.size(); i++)
    {
        auto[X_2,  Y_2]  = chi[1]->evaluate(sign2);
        auto[X_4,  Y_4]  = chi[3]->evaluate(sign4);
        auto[Xs_2, Ys_2] = chi[1]->derivative(sign2);
        auto[Xs_4, Ys_4] = chi[3]->derivative(sign4);
        Xs_2 *= sign2;
        Ys_2 *= sign2;
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        double x1i = x1(i);
        auto [X_1,  Y_1]  = chi[0]->evaluate(sign1*x1i);
        auto [X_3,  Y_3]  = chi[2]->evaluate(sign3*x1i);
        auto [Xs_3, Ys_3] = chi[2]->derivative(sign3*x1i);
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        dxdx1_north(i) = (X_2 - X_4 + 2*Xs_3)/2 - (x_3 - x_4)/2;
        dydx1_north(i) = (Y_2 - Y_4 + 2*Ys_3)/2 - (y_3 - y_4)/2;
        dxdx2_north(i) = ((1 - x1i)*Xs_4 + (1 + x1i)*Xs_2 + X_3 - X_1)/2 - ((1 - x1i)*(x_4 - x_1) + (1 + x1i)*(x_3 - x_2))/4;
        dydx2_north(i) = ((1 - x1i)*Ys_4 + (1 + x1i)*Ys_2 + Y_3 - Y_1)/2 - ((1 - x1i)*(y_4 - y_1) + (1 + x1i)*(y_3 - y_2))/4;
    }

    arma::vec    a_11_north = pow(dxdx1_north, 2) + pow(dydx1_north, 2) + pow(dzdx1.col(dzdx1.n_cols-1), 2);
    arma::vec    a_12_north = dxdx1_north%dxdx2_north + dydx1_north%dydx2_north + dzdx1.col(dzdx1.n_cols-1)%dzdx2.col(dzdx2.n_cols-1);
    arma::vec    a_22_north = pow(dxdx2_north, 2) + pow(dydx2_north, 2) + pow(dzdx2.col(dzdx2.n_cols-1), 2);

    arma::vec    a_north = a_11_north%a_22_north - pow(a_12_north, 2);

    arma::vec    a12_north =-a_12_north/a_north;
    arma::vec    a22_north = a_11_north/a_north;

    arma::vec    h_2s2_north = sqrt(a22_north);
    arma::vec    h_2s1_north = a12_north/sqrt(a22_north);

    return std::make_tuple(h_2s2_north, h_2s1_north);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsWest(const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactorsWest(x2, chi);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsWest(const arma::vec &x2, const std::array<Lagrange::CurveInterpolant*, 4> chi)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsWest", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);

    arma::rowvec dxdx1_west(x2.size(), arma::fill::none), dydx1_west(x2.size(), arma::fill::none),
                 dxdx2_west(x2.size(), arma::fill::none), dydx2_west(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto [X_1,  Y_1]  = chi[0]->evaluate(-sign1);
        auto [X_3,  Y_3]  = chi[2]->evaluate(-sign3);
        auto [Xs_1, Ys_1] = chi[0]->derivative(-sign1);
        auto [Xs_3, Ys_3] = chi[2]->derivative(-sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_4, Ys_4] = chi[3]->derivative(sign4*x2j);
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        dxdx1_west(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_west(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_west(j) = (2*Xs_4 + X_3 - X_1)/2 - (x_4 - x_1)/2;
        dydx2_west(j) = (2*Ys_4 + Y_3 - Y_1)/2 - (y_4 - y_1)/2;
    }

    arma::rowvec e_11_west = pow(dxdx1_west, 2) + pow(dydx1_west, 2);
    arma::rowvec e_12_west = dxdx1_west%dxdx2_west + dydx1_west%dydx2_west;
    arma::rowvec e_22_west = pow(dxdx2_west, 2) + pow(dydx2_west, 2);

    arma::rowvec e_west  = e_11_west %e_22_west  - pow(e_12_west,  2);

    arma::rowvec e11_west = e_22_west/e_west;
    arma::rowvec e12_west =-e_12_west/e_west;

    arma::rowvec h_1s1_west = sqrt(e11_west);
    arma::rowvec h_1s2_west = e12_west/sqrt(e11_west);

    return std::make_tuple(h_1s1_west, h_1s2_west);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsWest(const std::array<Lagrange::CurveInterpolant*, 4> chi, const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    arma::vec x2 = Chebyshev::gaussLobatto(chi[1]->getNodes().size());
    return covariantScaleFactorsWest(x2, chi, dzdx1, dzdx2);
}
std::tuple<arma::rowvec, arma::rowvec> Lagrange::covariantScaleFactorsWest(const arma::vec &x2, const std::array<Lagrange::CurveInterpolant*, 4> chi,
    const arma::mat &dzdx1, const arma::mat &dzdx2)
{
    auto [sign1, sign2, sign3, sign4] = signs("covariantScaleFactorsWest", chi);
    auto [x_1, y_1] = chi[0]->evaluate(-sign1);
    auto [x_2, y_2] = chi[0]->evaluate( sign1);
    auto [x_3, y_3] = chi[2]->evaluate( sign3);
    auto [x_4, y_4] = chi[2]->evaluate(-sign3);
    arma::rowvec dxdx1_west(x2.size(), arma::fill::none), dydx1_west(x2.size(), arma::fill::none),
                 dxdx2_west(x2.size(), arma::fill::none), dydx2_west(x2.size(), arma::fill::none);
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < x2.size(); j++)
    {
        auto [X_1,  Y_1]  = chi[0]->evaluate(-sign1);
        auto [X_3,  Y_3]  = chi[2]->evaluate(-sign3);
        auto [Xs_1, Ys_1] = chi[0]->derivative(-sign1);
        auto [Xs_3, Ys_3] = chi[2]->derivative(-sign3);
        Xs_1 *= sign1;
        Ys_1 *= sign1;
        Xs_3 *= sign3;
        Ys_3 *= sign3;
        double x2j = x2(j);
        auto [X_2,  Y_2]  = chi[1]->evaluate(sign2*x2j);
        auto [X_4,  Y_4]  = chi[3]->evaluate(sign4*x2j);
        auto [Xs_4, Ys_4] = chi[3]->derivative(sign4*x2j);
        Xs_4 *= sign4;
        Ys_4 *= sign4;
        dxdx1_west(j) = (X_2 - X_4 + (1 - x2j)*Xs_1 + (1 + x2j)*Xs_3)/2 - ((1 - x2j)*(x_2 - x_1) + (1 + x2j)*(x_3 - x_4))/4;
        dydx1_west(j) = (Y_2 - Y_4 + (1 - x2j)*Ys_1 + (1 + x2j)*Ys_3)/2 - ((1 - x2j)*(y_2 - y_1) + (1 + x2j)*(y_3 - y_4))/4;
        dxdx2_west(j) = (2*Xs_4 + X_3 - X_1)/2 - (x_4 - x_1)/2;
        dydx2_west(j) = (2*Ys_4 + Y_3 - Y_1)/2 - (y_4 - y_1)/2;
    }

    arma::rowvec a_11_west = pow(dxdx1_west, 2) + pow(dydx1_west, 2) + pow(dzdx1.row(0), 2);
    arma::rowvec a_12_west = dxdx1_west%dxdx2_west + dydx1_west%dydx2_west + dzdx1.row(0)%dzdx2.row(0);
    arma::rowvec a_22_west = pow(dxdx2_west, 2) + pow(dydx2_west, 2) + pow(dzdx2.row(0), 2);

    arma::rowvec a_west  = a_11_west %a_22_west  - pow(a_12_west,  2);

    arma::rowvec a11_west = a_22_west/a_west;
    arma::rowvec a12_west =-a_12_west/a_west;

    arma::rowvec h_1s1_west = sqrt(a11_west);
    arma::rowvec h_1s2_west = a12_west/sqrt(a11_west);

    return std::make_tuple(h_1s1_west, h_1s2_west);
}