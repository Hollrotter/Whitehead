#include <boost/math/special_functions/jacobi.hpp>
#include "JacobiAlpha.hpp"
#include "JacobiBeta.hpp"

int main()
{
    size_t n = 10;
    size_t m = 10;
    JacobiAlpha phiA(n, 0);
    JacobiBeta  phiB(m, 0);
    switch (1)
    {
        case 0: // Left and Right
        {
            double  TL   = phiA.constant();
            double  TR   = phiA.constant();
            double  TLp1 = phiA.linear(-1);
            double  TRp1 = phiA.linear( 1);
            double dTL   = phiA.constantDerivative();
            double dTR   = phiA.constantDerivative();
            double dTLp1 = phiA.linearDerivative();
            double dTRp1 = phiA.linearDerivative();
            for (size_t p = 0; p < n; p++)
            {
                double  tL = phiA.left(p);
                double  tR = phiA.right(p);
                double dtL = phiA.leftDerivative(p);
                double dtR = phiA.rightDerivative(p);

                std::cout << "  p = " <<   p << '\n';
                std::cout << " tL = " <<  tL << '\n';
                std::cout << " tR = " <<  tR << '\n';
                std::cout << "dtL = " << dtL << '\n';
                std::cout << "dtR = " << dtR << '\n';
                std::cout << "tL-TL = " << tL-TL << '\n';
                std::cout << "tR-TR = " << tR-TR << '\n';
                std::cout << "dtL-dTL = " << dtL-dTL << '\n';
                std::cout << "dtR-dTR = " << dtR-dTR << '\n';
                std::swap(TL, TLp1);
                std::swap(TR, TRp1);
                phiA.next(p, -1, TL, TLp1);
                phiA.next(p,  1, TR, TRp1);
                std::swap(dTL, dTLp1);
                std::swap(dTR, dTRp1);
                phiA.nextDerivative(p, -1, TL, dTL, dTLp1);
                phiA.nextDerivative(p,  1, TR, dTR, dTRp1);
            }
             TL   = phiB.constant();
             TR   = phiB.constant();
             TLp1 = phiB.linear(-1);
             TRp1 = phiB.linear( 1);
            dTL   = phiB.constantDerivative();
            dTR   = phiB.constantDerivative();
            dTLp1 = phiB.linearDerivative();
            dTRp1 = phiB.linearDerivative();
            for (size_t q = 0; q < m; q++)
            {
                double  tL = phiB.left(q);
                double  tR = phiB.right(q);
                double dtL = phiB.leftDerivative(q);
                double dtR = phiB.rightDerivative(q);

                std::cout << "  q = " <<   q << '\n';
                std::cout << " tL = " <<  tL << '\n';
                std::cout << " tR = " <<  tR << '\n';
                std::cout << "dtL = " << dtL << '\n';
                std::cout << "dtR = " << dtR << '\n';
                std::cout << "tL-TL = " << tL-TL << '\n';
                std::cout << "tR-TR = " << tR-TR << '\n';
                std::cout << "dtL-dTL = " << dtL-dTL << '\n';
                std::cout << "dtR-dTR = " << dtR-dTR << '\n';
                std::swap(TL, TLp1);
                std::swap(TR, TRp1);
                phiB.next(q, -1, TL, TLp1);
                phiB.next(q,  1, TR, TRp1);
                std::swap(dTL, dTLp1);
                std::swap(dTR, dTRp1);
                phiB.nextDerivative(q, -1, TL, dTL, dTLp1);
                phiB.nextDerivative(q,  1, TR, dTR, dTRp1);
            }
            break;
        }
        case 1:
        {
            size_t N = 5*n, M = 5*m;

            arma::mat t1(N, n);
            arma::mat t2(M, m);

            arma::vec x1 = arma::linspace(-1, 1, N);
            arma::vec x2 = arma::linspace(-1, 1, M);
            t1.col(0) = phiA.constant(N);
            t1.col(1) = phiA.linear(x1);
            for (size_t p = 1; p < n-1; p++)
                phiA.next(p, x1, t1);
            t2.col(0) = phiB.constant(N);
            t2.col(1) = phiB.linear(x2);
            for (size_t q = 1; q < m-1; q++)
                phiB.next(q, x2, t2);

            std::ofstream fileA("plot/Data/Jacobi/JacobiA");
            for (size_t i = 0; i < N; i++)
            {
                fileA << x1(i) << ' ';
                for (size_t p = 0; p < n; p++)
                    fileA << t1(i, p) << ' ';
                fileA << '\n';
            }
            fileA.close();
            std::ofstream fileB("plot/Data/Jacobi/JacobiB");
            for (size_t j = 0; j < M; j++)
            {
                fileB << x2(j) << ' ';
                for (size_t q = 0; q < m; q++)
                    fileB << t2(j, q) << ' ';
                fileB << '\n';
            }
            fileB.close();
            break;
        }
    }
}