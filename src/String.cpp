#include "String.hpp"

/**
 * @brief Construct a new String:: String object
 * 
 * @param _c 
 * @param _n 
 * @param _sigma 
 */
String::String(double _c, size_t _n, double _sigma) : c(_c), n(_n), sigma(_sigma)
{
    if (c <= 0 || n <= 0 || sigma <= 0)
    {
        std::println("All given parameters must be positive!");
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief 
 * 
 * @param _Et 
 */
void String::youngsModulus(double _Et)
{
    if (_Et < 0)
    {
        std::println("Young's modulus times thickness must be positive!");
        exit(EXIT_FAILURE);
    }
    Et = _Et;
}

/**
 * @brief 
 * 
 * @param _omega 
 */
void String::relaxationFactor(double _omega)
{
    try
    {
        if (_omega <= 0 || _omega > 1)
            throw std::runtime_error("Under-relaxation factor must lie between 0 and 1!");
    }
    catch (std::exception const &e)
    {
        std::cerr << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
    omega = _omega;
}

/**
 * @brief 
 * 
 * @param f 
 */
void String::load(std::function<double(double)> f)
{
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++)
        p(i) = f(x(i));
}

/**
 * @brief 
 * 
 * @param loc 
 * @param bc 
 * @param val 
 */
void String::boundary(Location loc, BC bc, double val)
{
    size_t i;
    int sign;
    switch(loc)
    {
        case Location::Front:
            i    = 0;
            sign = 1;
            front = bc;
            frontBC = val;
            break;
        case Location::Back:
            i    = n-1;
            sign =-1;
            back = bc;
            backBC = val;
            break;
    }
}

/**
 * @brief 
 * 
 * @param f 
 * @return double 
 */
double String::integrate(arma::vec f)
{
    f(0) = 0;
    arma::mat Dx = D1;
    Dx.row(0).eye();
    arma::vec int_f_dx = solve(Dx, f);
    return int_f_dx(n-1);
}

/**
 * @brief 
 * 
 * @param filename 
 */
void String::output(std::string filename)
{
    std::ofstream file(filename);
    for (size_t i = 0; i < n; i++)
        file << x(i) << ' ' << z(i) << '\n';
    file.close();
}