#include "Lagrange.hpp"

arma::vec Lagrange::CurveInterpolant::parametrize()
{
    return parametrize(Chebyshev::gaussLobatto(x.size()));
}

arma::vec Lagrange::CurveInterpolant::parametrize(arma::vec s)
{
    size_t n = x.size();
    arma::vec w = barycentricWeights(n);
    arma::vec t(n, arma::fill::zeros);
    double L;
    if (almostEqual(s(0),-1))
    {
        arma::mat D = Lagrange::derivativeMatrix(s);
        arma::vec dxds = D*x;
        arma::vec dyds = D*y;

        arma::vec I = sqrt(pow(dxds, 2) + pow(dyds, 2));
        D.row(0).eye();
        I(0) = 0;
        arma::vec int_I_dx = solve(D, I);
        L = int_I_dx(n-1);
    }
    else
    {
        for (size_t i = 0; i < n; i++)
        {
            fastgl::QuadPair gl = fastgl::GLPair(n, i+1);
            double dxds = interpolantDerivative(gl.x(), s, x, w);
            double dyds = interpolantDerivative(gl.x(), s, y, w);
            L += gl.weight * sqrt(pow(dxds, 2) + pow(dyds, 2));
        }
    }
    if (almostEqual(s(0),-1))
        t(0) =-1;
    else
    {
        double I1;
        double I2 = sqrt(pow(interpolantDerivative(-1, s, x, w), 2)
                       + pow(interpolantDerivative(-1, s, y, w), 2));
        t(0) = s(0);
        for (size_t q = 0; q < 1000; q++)
        {
            double dxdt = interpolantDerivative((t(0)-1)/2, s, x, w);
            double dydt = interpolantDerivative((t(0)-1)/2, s, y, w);

            double Im = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            dxdt = interpolantDerivative(t(0), s, x, w);
            dydt = interpolantDerivative(t(0), s, y, w);

            I1 = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            double f = L/2*(s(0) + 1) - (t(0) + 1)*(I1 + 4*Im + I2)/6;
            t(0) += 0.1*f/I1;
            if (fabs(f) < 1e-10)
                break;
        }
        I2 = I1;
    }
    double I2 = sqrt(pow(interpolantDerivative(-1, s, x, w), 2)
                   + pow(interpolantDerivative(-1, s, y, w), 2));
    double I1;
    for (size_t i = 1; i < n-1; i++)
    {
        t(i) = t(i-1) + s(i) - s(i-1);
        for (size_t q = 0; q < 1000; q++)
        {
            double dxdt = interpolantDerivative((t(i)+t(i-1))/2, s, x, w);
            double dydt = interpolantDerivative((t(i)+t(i-1))/2, s, y, w);

            double Im = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            dxdt = interpolantDerivative(t(i), s, x, w);
            dydt = interpolantDerivative(t(i), s, y, w);

            I1 = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            double f = L/2*(s(i) - s(i-1)) - (t(i) - t(i-1))*(I1 + 4*Im + I2)/6;
            t(i) += 0.1*f/I1;
            if (fabs(f) < 1e-10)
                break;
        }
        I2 = I1;
    }
    if (almostEqual(s(n-1), 1))
        t(n-1) = 1;
    else
    {
        t(n-1) = t(n-2) + s(n-1) - s(n-2);
        for (size_t q = 0; q < 1000; q++)
        {
            double dxdt = interpolantDerivative((t(n-1)+t(n-2))/2, s, x, w);
            double dydt = interpolantDerivative((t(n-1)+t(n-2))/2, s, y, w);

            double Im = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            dxdt = interpolantDerivative(t(n-1), s, x, w);
            dydt = interpolantDerivative(t(n-1), s, y, w);

            I1 = sqrt(pow(dxdt, 2) + pow(dydt, 2));

            double f = L/2*(s(n-1) - s(n-2)) - (t(n-1) - t(n-2))*(I1 + 4*Im + I2)/6;
            t(n-1) += 0.1*f/I1;
            if (fabs(f) < 1e-10)
                break;
        }
        I2 = I1;
    }
    return t;
}