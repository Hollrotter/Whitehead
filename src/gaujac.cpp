#include "gaujac.hpp"

std::pair<arma::vec, arma::vec> gaujac(size_t n, double alpha, double beta)
{
    arma::vec x(n), w(n);
    const size_t MAXIT = 10;
    const double EPS = 1.0e-14; // EPS is the relative precision.
    double z = 0;
    for (size_t i = 0; i < n; i++) // Loop over the desired roots.
    {
        if (i == 0) // Initial guess for the largest root.
        {
            double an = alpha/n;
            double bn = beta/n;
            double r1 = (1+alpha)*(2.78/(4+n*n)+0.768*an/n);
            double r2 = 1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
            z = 1-r1/r2;
        }
        else if (i == 1) // Initial guess for the second largest root.
        {
            double r1 = (4.1+alpha)/((1+alpha)*(1+0.156*alpha));
            double r2 = 1+0.06*(n-8)*(1+0.12*alpha)/n;
            double r3 = 1+0.012*beta*(1+0.25*fabs(alpha))/n;
            z -= (1-z)*r1*r2*r3;
        }
        else if (i == 2) // Initial guess for the third largest root.
        {
            double r1 = (1.67+0.28*alpha)/(1+0.37*alpha);
            double r2 = 1+0.22*(1-8./n);
            double r3 = 1+8*beta/((6.28+beta)*n*n);
            z -= (x(0)-z)*r1*r2*r3;
        }
        else if (i == n-2) // Initial guess for the second smallest root.
        {
            double r1 = (1+0.235*beta)/(0.766+0.119*beta);
            double r2 = 1/(1+0.639*(n-4)/(1+0.71*(n-4)));
            double r3 = 1/(1+20*alpha/((7.5+alpha)*n*n));
            z += (z-x(n-4))*r1*r2*r3;
        }
        else if (i == n-1) // Initial guess for the smallest root.
        {
            double r1 = (1+0.37*beta)/(1.67+0.28*beta);
            double r2 = 1/(1+0.22*(n-8)/n);
            double r3 = 1/(1+8*alpha/((6.28+alpha)*n*n));
            z += (z-x(n-3))*r1*r2*r3;
        }
        else // Initial guess for the other roots.
            z = 3*x(i-1)-3*x(i-2)+x(i-3);
        double alphabeta = alpha+beta;
        double temp, p2, pp;
        for (size_t its = 1; its <= MAXIT; its++) // Refinement by Newton's method.
        {
            // Start the recurrence with P0 and P1 to avoid a division by zero when alpha+beta=0 or -1.
            temp = 2 + alphabeta; 
            double p1 = (alpha-beta+temp*z)/2;
            p2 = 1;
            // Loop up the recurrence relation to get the Jacobi polynomial evaluated at z.
            for (size_t j = 2; j <= n; j++)
            {
                double p3 = p2;
                p2 = p1;
                temp = 2*j+alphabeta;
                double a = 2*j*(j+alphabeta)*(temp-2);
                double b = (temp-1)*(alpha*alpha-beta*beta+temp*(temp-2)*z);
                double c = 2*(j-1+alpha)*(j-1+beta)*temp;
                p1 = (b*p2-c*p3)/a;
            }
            pp = (n*(alpha-beta-temp*z)*p1+2*(n+alpha)*(n+beta)*p2)/(temp*(1-z*z));
            // p1 is now the desired Jacobi polynomial. We next compute pp, its derivative, by
            // a standard relation involving also p2, the polynomial of one lower order.
            double z1 = z;
            z = z1 - p1/pp; // Newton's formula.
            if (fabs(z-z1) <= EPS)
                break;
            if (its == MAXIT)
                throw("Too many iterations in gaujac!");
        }
        x(i) = z;
        w(i) = exp(lgamma(alpha+n) + lgamma(beta+n) - lgamma(n+1) - lgamma(n+alphabeta+1))
             * temp*pow(2, alphabeta)/(pp*p2);
    }
    return {reverse(x), reverse(w)};
}