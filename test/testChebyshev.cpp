#include "Chebyshev.hpp"

int main()
{
    arma::vec x   = Chebyshev::gaussLobatto(100);
    arma::vec z   = exp(-pow(3*x, 2));
    arma::vec dz  = Chebyshev::derivativeMatrix(x, Derivative::first)  * z;
    arma::vec d2z = Chebyshev::derivativeMatrix(x, Derivative::second) * z;
    z.t().print();
    dz.t().print();
    d2z.t().print();
    std::ofstream file("plot/Data/Chebyshev");
    for (size_t i = 0; i < x.size(); i++)
        file << x(i) << ' ' << z(i) << ' ' << dz(i) << ' ' << d2z(i) << '\n';
    file.close();
}