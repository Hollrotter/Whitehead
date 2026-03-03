#include "B_Spline.hpp"

/**
 * @brief Construct a new b spline::b spline object
 * 
 * @param v 
 * @param _order Order of the Spline.
 */
B_Spline::B_Spline(arma::vec v, size_t _order) : order(_order), X(arma::zeros(v.size()+2*order))
{
    if (order == 0)
        X = v;
    else
    {
        X.head(order).fill(v(0));
        X.subvec(order, X.size()-1-order) = v;
        X.tail(order).fill(v(v.size()-1));
    }
}

arma::vec B_Spline::operator()(const double x)
{
    arma::vec y(X.size()-1-order);
    #pragma omp parallel for
    for (size_t j = 0; j < y.size(); j++)
        y(j) = operator()(x, j, order);
    return y;
}

/**
 * @brief 
 * 
 * @param x 
 * @param i 
 * @param k 
 * @return double 
 */
double B_Spline::operator()(const double x, const size_t i, const size_t k)
{
    if (k == 0)
        if (almostEqual(x, X.back()))
            return (x > X(i) || almostEqual(x, X(i))) && (x < X(i+1) || almostEqual(x, X(i+1))) ? 1 : 0;
        else
            return (x > X(i) || almostEqual(x, X(i))) && (x < X(i+1)) ? 1 : 0;
    double c1 = (almostEqual(X(i+k),   X(i))   == false) ? (x - X(i))    /(X(k+i) - X(i))    *operator()(x, i,   k-1) : 0;
    double c2 = (almostEqual(X(i+k+1), X(i+1)) == false) ? (x - X(i+k+1))/(X(i+1) - X(i+k+1))*operator()(x, i+1, k-1) : 0;
    return c1 + c2;
}

arma::mat B_Spline::operator()(const arma::vec x)
{
    arma::mat y(X.size()-1-order, x.size());
    for (size_t i = 0; i < x.size(); i++)
        y.col(i) = operator()(x(i));
    return y;
}

/**
 * @brief 
 * 
 * @param x 
 * @param i 
 * @param k 
 * @return arma::vec 
 */
arma::vec B_Spline::operator()(const arma::vec x, const size_t i, const size_t k)
{
    arma::vec y(x.size());
    #pragma omp parallel for
    for (size_t j = 0; j < x.size(); j++)
        y(j) = operator()(x(j), i, k);
    return y;
}

arma::vec B_Spline::diff(const double x)
{
    arma::vec dy(X.size()-1-order);
    #pragma omp parallel for
    for (size_t j = 0; j < dy.size(); j++)
        dy(j) = diff(x, j, order);
    return dy;
}

/**
 * @brief 
 * 
 * @param x 
 * @param i 
 * @param k 
 * @return double 
 */
double B_Spline::diff(const double x, const size_t i, const size_t k)
{
    double c1 = (almostEqual(X(i+k),   X(i))   == false) ? operator()(x, i,   k-1)/(X(i+k)   - X(i))   : 0;
    double c2 = (almostEqual(X(i+k+1), X(i+1)) == false) ? operator()(x, i+1, k-1)/(X(i+k+1) - X(i+1)) : 0;
    return k*(c1 - c2);
}

/**
 * @brief 
 * 
 * @param x 
 * @param i 
 * @param k 
 * @return arma::vec 
 */
arma::vec B_Spline::diff(const arma::vec x, const size_t i, const size_t k)
{
    arma::vec y(x.size());
    #pragma omp parallel for
    for (size_t j = 0; j < x.size(); j++)
        y(j) = diff(x(j), i, k);
    return y;
}

double B_Spline::integrate(const size_t i, const size_t j)
{
    double integral = 0;
    B_Spline B(X.subvec(order, X.size()-1-order), order+1);
    #pragma omp parallel for reduction (+:integral)
    for (size_t n = j+1; n < B.X.size()-2*B.order+order; n++)
        integral += B(B.X(i+1+B.order), n, order+1) - B(B.X(i+B.order), n, order+1);
    return (X(j+order+1) - X(j))/(order+1)*integral;
}

double B_Spline::integrate(const size_t i, const size_t j, const size_t k)
{
    double integral = 0;
    #pragma omp parallel for reduction (+:integral)
    for (size_t n = j; n < X.size()-order-1; n++)
        integral += operator()(X(i+1+order), n, k+1) - operator()(X(i+order), n, k+1);
    return (X(j+k+1) - X(j))/(k+1)*integral;
}

double B_Spline::integrate(const double x1, const double x2, const size_t j)
{
    double integral = 0;
    B_Spline B(X.subvec(order, X.size()-1-order), order+1);
    #pragma omp parallel for reduction (+:integral)
    for (size_t n = j+1; n < B.X.size()-2*B.order+order; n++)
        integral += B(x2, n, order+1) - B(x1, n, order+1);
    return (X(j+order+1) - X(j))/(order+1)*integral;
}

double B_Spline::integrate(const double x1, const double x2, const size_t j, const size_t k)
{
    double integral = 0;
    #pragma omp parallel for reduction (+:integral)
    for (size_t n = j; n < X.size()-order-1; n++)
        integral += operator()(x2, n, k+1) - operator()(x1, n, k+1);
    return (X(j+k+1) - X(j))/(k+1)*integral;
}

double B_Spline::integrate_x(const double x1, const double x2, const size_t j)
{
    double integral = 0;
    B_Spline B1(X.subvec(order, X.size()-1-order), order+1);
    #pragma omp parallel for reduction (+:integral)
    for (size_t n = j+1; n < B1.X.size()-2*B1.order+order; n++)
        integral += x2*B1(x2, n) - x1*B1(x1, n) - B1.integrate(x1, x2, n);
    return (X(j+order+1) - X(j))/(order+1)*integral;
}