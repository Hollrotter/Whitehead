#include "Wing.hpp"

void Wing::integralRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xp)* x1 + xp+1)/2;
            dx1_gl = ((1-xp)*dx1 + xp+1)/2;
            w1_gl  = pow((1-xp)/2, 1.5)* w1%sqrt(1+x1_gl);
            dw1_gl = pow((1-xp)/2, 0.5)*dw1/sqrt(1+x1_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xp)* x1 + xp+1)/2;
            dx1_gl = ((1-xp)*dx1 + xp+1)/2;
            w1_gl  = pow((1-xp)/2, 1.5)* w1;
            dw1_gl = pow((1-xp)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (1 - xp)/2 + xp;
                w1_gl(ii) = gl_x.weight*(1 - xp)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    arma::vec  x2_gl = phi2->xg;
    arma::vec  w2_gl = phi2->wg;
    arma::vec dx2_gl = (phi2->basis == Basis::T) ? x2_gl : phi2->dxg;
    arma::vec dw2_gl = phi2->dwg;

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xm+1)* x1 + xm-1)/2;
            dx1_gl = ((xm+1)*dx1 + xm-1)/2;
            w1_gl  = pow((xm+1)/2, 1.5)* w1%sqrt(1-x1_gl);
            dw1_gl = pow((xm+1)/2, 0.5)*dw1/sqrt(1-x1_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xm+1)* x1 + xm-1)/2;
            dx1_gl = ((xm+1)*dx1 + xm-1)/2;
            w1_gl  = pow((xm+1)/2, 1.5)* w1;
            dw1_gl = pow((xm+1)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (xm + 1)/2 - 1;
                w1_gl(ii) = gl_x.weight*(xm + 1)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    arma::vec  x2_gl = phi2->xg;
    arma::vec  w2_gl = phi2->wg;
    arma::vec dx2_gl = (phi2->basis == Basis::T) ? x2_gl : phi2->dxg;
    arma::vec dw2_gl = phi2->dwg;

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralUpperNonlinear(size_t i, size_t j, size_t n1, size_t n2, double yp)
{
    arma::vec  x1_gl = phi1->xg;
    arma::vec  w1_gl = phi1->wg;
    arma::vec dx1_gl = (phi1->basis == Basis::T) ? x1_gl : phi1->dxg;
    arma::vec dw1_gl = phi1->dwg;
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-yp)* x2 + yp+1)/2;
            dx2_gl = ((1-yp)*dx2 + yp+1)/2;
            w2_gl  = pow((1-yp)/2, 1.5)* w2%sqrt(1+x2_gl);
            dw2_gl = pow((1-yp)/2, 0.5)*dw2/sqrt(1+x2_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-yp)* x2 + yp+1)/2;
            dx2_gl = ((1-yp)*dx2 + yp+1)/2;
            w2_gl  = pow((1-yp)/2, 1.5)* w2;
            dw2_gl = pow((1-yp)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (1 - yp)/2 + yp;
                w2_gl(jj) = gl_y.weight*(1 - yp)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralLowerNonlinear(size_t i, size_t j, size_t n1, size_t n2, double ym)
{
    arma::vec  x1_gl = phi1->xg;
    arma::vec  w1_gl = phi1->wg;
    arma::vec dx1_gl = (phi1->basis == Basis::T) ? x1_gl : phi1->dxg;
    arma::vec dw1_gl = phi1->dwg;
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((ym+1)* x2 + ym-1)/2;
            dx2_gl = ((ym+1)*dx2 + ym-1)/2;
            w2_gl  = pow((ym+1)/2, 1.5)* w2%sqrt(1-x2_gl);
            dw2_gl = pow((ym+1)/2, 0.5)*dw2/sqrt(1-x2_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((ym+1)* x2 + ym-1)/2;
            dx2_gl = ((ym+1)*dx2 + ym-1)/2;
            w2_gl  = pow((ym+1)/2, 1.5)* w2;
            dw2_gl = pow((ym+1)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (ym + 1)/2 - 1;
                w2_gl(jj) = gl_y.weight*(ym + 1)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralLowerRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp, double yp)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xp)* x1 + xp+1)/2;
            dx1_gl = ((1-xp)*dx1 + xp+1)/2;
            w1_gl  = pow((1-xp)/2, 1.5)* w1%sqrt(1+x1_gl);
            dw1_gl = pow((1-xp)/2, 0.5)*dw1/sqrt(1+x1_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xp)* x1 + xp+1)/2;
            dx1_gl = ((1-xp)*dx1 + xp+1)/2;
            w1_gl  = pow((1-xp)/2, 1.5)* w1;
            dw1_gl = pow((1-xp)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (1 - xp)/2 + xp;
                w1_gl(ii) = gl_x.weight*(1 - xp)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((yp+1)* x2 + yp-1)/2;
            dx2_gl = ((yp+1)*dx2 + yp-1)/2;
            w2_gl  = pow((yp+1)/2, 1.5)* w2%sqrt(1-x2_gl);
            dw2_gl = pow((yp+1)/2, 0.5)*dw2/sqrt(1-x2_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((yp+1)* x2 + yp-1)/2;
            dx2_gl = ((yp+1)*dx2 + yp-1)/2;
            w2_gl  = pow((yp+1)/2, 1.5)* w2;
            dw2_gl = pow((yp+1)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (yp + 1)/2 - 1;
                w2_gl(jj) = gl_y.weight*(yp + 1)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralLowerLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xp, double ym)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xp+1)* x1 + xp-1)/2;
            dx1_gl = ((xp+1)*dx1 + xp-1)/2;
            w1_gl  = pow((xp+1)/2, 1.5)* w1%sqrt(1-x1_gl);
            dw1_gl = pow((xp+1)/2, 0.5)*dw1/sqrt(1-x1_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xp+1)* x1 + xp-1)/2;
            dx1_gl = ((xp+1)*dx1 + xp-1)/2;
            w1_gl  = pow((xp+1)/2, 1.5)* w1;
            dw1_gl = pow((xp+1)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (xp + 1)/2 - 1;
                w1_gl(ii) = gl_x.weight*(xp + 1)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((ym+1)* x2 + ym-1)/2;
            dx2_gl = ((ym+1)*dx2 + ym-1)/2;
            w2_gl  = pow((ym+1)/2, 1.5)* w2%sqrt(1-x2_gl);
            dw2_gl = pow((ym+1)/2, 0.5)*dw2/sqrt(1-x2_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACB[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACB2[n2]);
            x2_gl  = ((ym+1)* x2 + ym-1)/2;
            dx2_gl = ((ym+1)*dx2 + ym-1)/2;
            w2_gl  = pow((ym+1)/2, 1.5)* w2;
            dw2_gl = pow((ym+1)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (ym + 1)/2 - 1;
                w2_gl(jj) = gl_y.weight*(ym + 1)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralUpperRightNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double yp)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xm)* x1 + xm+1)/2;
            dx1_gl = ((1-xm)*dx1 + xm+1)/2;
            w1_gl  = pow((1-xm)/2, 1.5)* w1%sqrt(1+x1_gl);
            dw1_gl = pow((1-xm)/2, 0.5)*dw1/sqrt(1+x1_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACA[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACA2[n1]);
            x1_gl  = ((1-xm)* x1 + xm+1)/2;
            dx1_gl = ((1-xm)*dx1 + xm+1)/2;
            w1_gl  = pow((1-xm)/2, 1.5)* w1;
            dw1_gl = pow((1-xm)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (1 - xm)/2 + xm;
                w1_gl(ii) = gl_x.weight*(1 - xm)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-yp)* x2 + yp+1)/2;
            dx2_gl = ((1-yp)*dx2 + yp+1)/2;
            w2_gl  = pow((1-yp)/2, 1.5)* w2%sqrt(1+x2_gl);
            dw2_gl = pow((1-yp)/2, 0.5)*dw2/sqrt(1+x2_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-yp)* x2 + yp+1)/2;
            dx2_gl = ((1-yp)*dx2 + yp+1)/2;
            w2_gl  = pow((1-yp)/2, 1.5)* w2;
            dw2_gl = pow((1-yp)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (1 - yp)/2 + yp;
                w2_gl(jj) = gl_y.weight*(1 - yp)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralUpperLeftNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double ym)
{
    arma::vec x1_gl(n1, arma::fill::none), dx1_gl(n1, arma::fill::none);
    arma::vec x2_gl(n2, arma::fill::none), dx2_gl(n2, arma::fill::none);
    arma::vec w1_gl(n1, arma::fill::none), dw1_gl(n1, arma::fill::none);
    arma::vec w2_gl(n2, arma::fill::none), dw2_gl(n2, arma::fill::none);
    switch (phi1->basis)
    {
        case Basis::U:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xm+1)* x1 + xm-1)/2;
            dx1_gl = ((xm+1)*dx1 + xm-1)/2;
            w1_gl  = pow((xm+1)/2, 1.5)* w1%sqrt(1-x1_gl);
            dw1_gl = pow((xm+1)/2, 0.5)*dw1/sqrt(1-x1_gl);
            break;
        }
        case Basis::PB:
        {
            auto [ x1,  w1] = convert<NMAX>(n1, GAUJACB[n1]);
            auto [dx1, dw1] = convert<NMAX>(n1, GAUJACB2[n1]);
            x1_gl  = ((xm+1)* x1 + xm-1)/2;
            dx1_gl = ((xm+1)*dx1 + xm-1)/2;
            w1_gl  = pow((xm+1)/2, 1.5)* w1;
            dw1_gl = pow((xm+1)/2, 0.5)*dw1;
            break;
        }
        default:
        {
            for (size_t ii = 0; ii < n1; ii++)
            {
                fastgl::QuadPair gl_x = fastgl::GLPair(n1, ii+1);
                x1_gl(ii) = (1 - gl_x.x()) * (xm + 1)/2 - 1;
                w1_gl(ii) = gl_x.weight*(xm + 1)/2;
            }
            dx1_gl = x1_gl;
            if (phi1->basis == Basis::T)
                dw1_gl.zeros();
            else
            {
                dw1_gl = w1_gl;
                w1_gl  %= phi1->weightFunction(x1_gl);
                dw1_gl %= phi1->weightFunctionDerivative(x1_gl);
            }
            break;
        }
    }
    switch (phi2->basis)
    {
        case Basis::U:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-ym)* x2 + ym+1)/2;
            dx2_gl = ((1-ym)*dx2 + ym+1)/2;
            w2_gl  = pow((1-ym)/2, 1.5)* w2%sqrt(1+x2_gl);
            dw2_gl = pow((1-ym)/2, 0.5)*dw2/sqrt(1+x2_gl);
            break;
        }
        case Basis::PA:
        {
            auto [ x2,  w2] = convert<NMAX>(n2, GAUJACA[n2]);
            auto [dx2, dw2] = convert<NMAX>(n2, GAUJACA2[n2]);
            x2_gl  = ((1-ym)* x2 + ym+1)/2;
            dx2_gl = ((1-ym)*dx2 + ym+1)/2;
            w2_gl  = pow((1-ym)/2, 1.5)* w2;
            dw2_gl = pow((1-ym)/2, 0.5)*dw2;
            break;
        }
        default:
        {
            for (size_t jj = 0; jj < n2; jj++)
            {
                fastgl::QuadPair gl_y = fastgl::GLPair(n2, jj+1);
                x2_gl(jj) = (1 - gl_y.x()) * (1 - ym)/2 + ym;
                w2_gl(jj) = gl_y.weight*(1 - ym)/2;
            }
            dx2_gl = x2_gl;
            if (phi2->basis == Basis::T)
                dw2_gl.zeros();
            else
            {
                dw2_gl = w2_gl;
                w2_gl  %= phi2->weightFunction(x2_gl);
                dw2_gl %= phi2->weightFunctionDerivative(x2_gl);
            }
            break;
        }
    }

    integralNonlinearEval(i, j, n1, n2, x1_gl, x2_gl, dx1_gl, dx2_gl, w1_gl, w2_gl, dw1_gl, dw2_gl);
}

void Wing::integralNonlinearEval(size_t &i, size_t &j, size_t &n1, size_t &n2,
    arma::vec &x1_gl, arma::vec &x2_gl, arma::vec &dx1_gl, arma::vec &dx2_gl,
    arma::vec &w1_gl, arma::vec &w2_gl, arma::vec &dw1_gl, arma::vec &dw2_gl)
{
    arma::mat Tx_gl_0 = Lagrange::interpolationMatrix(x1,  x1_gl);
    arma::mat Ty_gl_0 = Lagrange::interpolationMatrix(x2,  x2_gl);
    arma::mat Tx_gl_1 = Lagrange::interpolationMatrix(x1, dx1_gl);
    arma::mat Ty_gl_2 = Lagrange::interpolationMatrix(x2, dx2_gl);
    arma::mat z_0_gl  = Lagrange::interpolation2D(Tx_gl_0, Ty_gl_0, z,  x1_gl,  x2_gl);
    arma::mat z_1_gl  = Lagrange::interpolation2D(Tx_gl_1, Ty_gl_0, z, dx1_gl,  x2_gl);
    arma::mat z_2_gl  = Lagrange::interpolation2D(Tx_gl_0, Ty_gl_2, z,  x1_gl, dx2_gl);
    arma::mat D1_gl_0 = Lagrange::derivativeMatrix(x1_gl);
    arma::mat D2_gl_0 = Lagrange::derivativeMatrix(x2_gl);
    arma::mat D1_gl_1 = Lagrange::derivativeMatrix(dx1_gl);
    arma::mat D2_gl_2 = Lagrange::derivativeMatrix(dx2_gl);
    arma::mat dz_0_gldx1 = D1_gl_0*z_0_gl;
    arma::mat dz_0_gldx2 = z_0_gl*D2_gl_0.t();
    arma::mat dz_1_gldx1 = D1_gl_1*z_1_gl;
    arma::mat dz_1_gldx2 = z_1_gl*D2_gl_0.t();
    arma::mat dz_2_gldx1 = D1_gl_0*z_2_gl;
    arma::mat dz_2_gldx2 = z_2_gl*D2_gl_2.t();
    auto [x_0_gl, y_0_gl] = Lagrange::TransfiniteQuadMap( x1_gl,  x2_gl, chi);
    auto [x_1_gl, y_1_gl] = Lagrange::TransfiniteQuadMap(dx1_gl,  x2_gl, chi);
    auto [x_2_gl, y_2_gl] = Lagrange::TransfiniteQuadMap( x1_gl, dx2_gl, chi);
    auto [dx_0_gldx1, dx_0_gldx2, dy_0_gldx1, dy_0_gldx2] = Lagrange::TransfiniteQuadMetrics( x1_gl,  x2_gl, chi);
    auto [dx_1_gldx1, dx_1_gldx2, dy_1_gldx1, dy_1_gldx2] = Lagrange::TransfiniteQuadMetrics(dx1_gl,  x2_gl, chi);
    auto [dx_2_gldx1, dx_2_gldx2, dy_2_gldx1, dy_2_gldx2] = Lagrange::TransfiniteQuadMetrics( x1_gl, dx2_gl, chi);
    arma::field<arma::mat> J0_gl = {{dx_0_gldx1, dx_0_gldx2}, {dy_0_gldx1, dy_0_gldx2}};
    arma::field<arma::mat> J1_gl = {{dx_1_gldx1, dx_1_gldx2}, {dy_1_gldx1, dy_1_gldx2}};
    arma::field<arma::mat> J2_gl = {{dx_2_gldx1, dx_2_gldx2}, {dy_2_gldx1, dy_2_gldx2}};
    arma::cube e_c_gl_0 = MetricCo(J0_gl);
    arma::cube e_c_gl_1 = MetricCo(J1_gl);
    arma::cube e_c_gl_2 = MetricCo(J2_gl);
    arma::cube ec_gl_0  = MetricContra(e_c_gl_0);
    arma::cube ec_gl_1  = MetricContra(e_c_gl_1);
    arma::cube ec_gl_2  = MetricContra(e_c_gl_2);
    arma::mat e_gl_0 = e_c_gl_0.slice(0)%e_c_gl_0.slice(2) - pow(e_c_gl_0.slice(1), 2);
    arma::mat e_gl_1 = e_c_gl_1.slice(0)%e_c_gl_1.slice(2) - pow(e_c_gl_1.slice(1), 2);
    arma::mat e_gl_2 = e_c_gl_2.slice(0)%e_c_gl_2.slice(2) - pow(e_c_gl_2.slice(1), 2);
    arma::mat sqrt_a0 = sqrt(e_gl_0%(1 + ec_gl_0.slice(0)%pow(dz_0_gldx1, 2) + 2*ec_gl_0.slice(1)%dz_0_gldx1%dz_0_gldx2 + ec_gl_0.slice(2)%pow(dz_0_gldx2, 2)));
    arma::mat sqrt_a1 = sqrt(e_gl_1%(1 + ec_gl_1.slice(0)%pow(dz_1_gldx1, 2) + 2*ec_gl_1.slice(1)%dz_1_gldx1%dz_1_gldx2 + ec_gl_1.slice(2)%pow(dz_1_gldx2, 2)));
    arma::mat sqrt_a2 = sqrt(e_gl_2%(1 + ec_gl_2.slice(0)%pow(dz_2_gldx1, 2) + 2*ec_gl_2.slice(1)%dz_2_gldx1%dz_2_gldx2 + ec_gl_2.slice(2)%pow(dz_2_gldx2, 2)));
    for (size_t ii = 0; ii < n1; ii++)
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r_0 = {x_0_gl(ii, jj) - xC(i, j), y_0_gl(ii, jj) - yC(i, j), z_0_gl(ii, jj) - zC(i, j)};
            arma::vec::fixed<3> r_1 = {x_1_gl(ii, jj) - xC(i, j), y_1_gl(ii, jj) - yC(i, j), z_1_gl(ii, jj) - zC(i, j)};
            arma::vec::fixed<3> r_2 = {x_2_gl(ii, jj) - xC(i, j), y_2_gl(ii, jj) - yC(i, j), z_2_gl(ii, jj) - zC(i, j)};
            double r3_0 = pow(norm(r_0), 3);
            double r3_1 = pow(norm(r_1), 3);
            double r3_2 = pow(norm(r_2), 3);
            arma::vec::fixed<3> n0_gl = arma::vec::fixed<3>({dy_0_gldx1(ii, jj)*dz_0_gldx2(ii, jj)-dz_0_gldx1(ii, jj)*dy_0_gldx2(ii, jj),
                                                             dz_0_gldx1(ii, jj)*dx_0_gldx2(ii, jj)-dx_0_gldx1(ii, jj)*dz_0_gldx2(ii, jj),
                                                             dx_0_gldx1(ii, jj)*dy_0_gldx2(ii, jj)-dy_0_gldx1(ii, jj)*dx_0_gldx2(ii, jj)})/sqrt_a0(ii, jj);
            arma::vec::fixed<3> n1_gl = arma::vec::fixed<3>({dy_1_gldx1(ii, jj)*dz_1_gldx2(ii, jj)-dz_1_gldx1(ii, jj)*dy_1_gldx2(ii, jj),
                                                             dz_1_gldx1(ii, jj)*dx_1_gldx2(ii, jj)-dx_1_gldx1(ii, jj)*dz_1_gldx2(ii, jj),
                                                             dx_1_gldx1(ii, jj)*dy_1_gldx2(ii, jj)-dy_1_gldx1(ii, jj)*dx_1_gldx2(ii, jj)})/sqrt_a1(ii, jj);
            arma::vec::fixed<3> n2_gl = arma::vec::fixed<3>({dy_2_gldx1(ii, jj)*dz_2_gldx2(ii, jj)-dz_2_gldx1(ii, jj)*dy_2_gldx2(ii, jj),
                                                             dz_2_gldx1(ii, jj)*dx_2_gldx2(ii, jj)-dx_2_gldx1(ii, jj)*dz_2_gldx2(ii, jj),
                                                             dx_2_gldx1(ii, jj)*dy_2_gldx2(ii, jj)-dy_2_gldx1(ii, jj)*dx_2_gldx2(ii, jj)})/sqrt_a2(ii, jj);
            arma::mat::fixed<3, 2> J0_red = {{dy_0_gldx2(ii, jj)*n0_gl(2) - n0_gl(1)*dz_0_gldx2(ii, jj),n0_gl(1)*dz_0_gldx1(ii, jj) - dy_0_gldx1(ii, jj)*n0_gl(2)},
                                             {dz_0_gldx2(ii, jj)*n0_gl(0) - n0_gl(2)*dx_0_gldx2(ii, jj),n0_gl(2)*dx_0_gldx1(ii, jj) - dz_0_gldx1(ii, jj)*n0_gl(0)},
                                             {dx_0_gldx2(ii, jj)*n0_gl(1) - n0_gl(0)*dy_0_gldx2(ii, jj),n0_gl(0)*dy_0_gldx1(ii, jj) - dx_0_gldx1(ii, jj)*n0_gl(1)}};
            arma::mat::fixed<3, 2> J1_red = {{dy_1_gldx2(ii, jj)*n1_gl(2) - n1_gl(1)*dz_1_gldx2(ii, jj),n1_gl(1)*dz_1_gldx1(ii, jj) - dy_1_gldx1(ii, jj)*n1_gl(2)},
                                             {dz_1_gldx2(ii, jj)*n1_gl(0) - n1_gl(2)*dx_1_gldx2(ii, jj),n1_gl(2)*dx_1_gldx1(ii, jj) - dz_1_gldx1(ii, jj)*n1_gl(0)},
                                             {dx_1_gldx2(ii, jj)*n1_gl(1) - n1_gl(0)*dy_1_gldx2(ii, jj),n1_gl(0)*dy_1_gldx1(ii, jj) - dx_1_gldx1(ii, jj)*n1_gl(1)}};
            arma::mat::fixed<3, 2> J2_red = {{dy_2_gldx2(ii, jj)*n2_gl(2) - n2_gl(1)*dz_2_gldx2(ii, jj),n2_gl(1)*dz_2_gldx1(ii, jj) - dy_2_gldx1(ii, jj)*n2_gl(2)},
                                             {dz_2_gldx2(ii, jj)*n2_gl(0) - n2_gl(2)*dx_2_gldx2(ii, jj),n2_gl(2)*dx_2_gldx1(ii, jj) - dz_2_gldx1(ii, jj)*n2_gl(0)},
                                             {dx_2_gldx2(ii, jj)*n2_gl(1) - n2_gl(0)*dy_2_gldx2(ii, jj),n2_gl(0)*dy_2_gldx1(ii, jj) - dx_2_gldx1(ii, jj)*n2_gl(1)}};
            double t2   = phi2->constant(); // First Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Polynomial (j+1 -> jp1)

            double t2_2   = phi2->constant(); // First Polynomial (j)
            double t2p1_2 = phi2->linear(dx2_gl(jj)); // Second Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Polynomial 2-direction
            {
                double t1   = phi1->constant(); // First Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Polynomial (i+1 -> ip1)

                double t1_1   = phi1->constant(); // First Polynomial (i)
                double t1p1_1 = phi1->linear(dx1_gl(ii)); // Second Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Polynomial 1-direction
                {
                    double dmu0dx1 = dt1 * t2;
                    double dmu0dx2 = t1 * dt2;
                    double dmu1dx1 = t1_1 * t2 * phi1->weightFunctionDerivativeFactor(dx1_gl(ii));
                    double dmu2dx2 = t1 * t2_2 * phi2->weightFunctionDerivativeFactor(dx2_gl(jj));
                    arma::vec::fixed<2> dmu0dxi = {dmu0dx1, dmu0dx2};
                    arma::vec::fixed<2> dmu1dxi = {dmu1dx1, dmu0dx2};
                    arma::vec::fixed<2> dmu2dxi = {dmu0dx1, dmu2dx2};
                    arma::vec::fixed<3> gradmu0 = J0_red*dmu0dxi;
                    arma::vec::fixed<3> gradmu1 = J1_red*dmu1dxi;
                    arma::vec::fixed<3> gradmu2 = J2_red*dmu2dxi;
                    arma::vec::fixed<3> gamma_gl_0 = cross(gradmu0, n0_gl);
                    arma::vec::fixed<3> gamma_gl_1 = cross(gradmu1, n1_gl);
                    arma::vec::fixed<3> gamma_gl_2 = cross(gradmu2, n2_gl);
                    arma::vec::fixed<3> q0_mu =  w1_gl(ii) *  w2_gl(jj) * cross(gamma_gl_0, r_0)/r3_0;
                    arma::vec::fixed<3> q1_mu = dw1_gl(ii) *  w2_gl(jj) * cross(gamma_gl_1, r_1)/r3_1;
                    arma::vec::fixed<3> q2_mu =  w1_gl(ii) * dw2_gl(jj) * cross(gamma_gl_2, r_2)/r3_2;

                    A(i+j*nx, p+q*nx) += dot(q0_mu+q1_mu+q2_mu, nC.row(i+j*nx));

                    std::swap(t1, t1p1); // Swap order so the oldest Polynomial will be overwritten
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(t1_1, t1p1_1); // Swap order so the oldest Polynomial will be overwritten
                    phi1->next(p, dx1_gl(ii), t1_1, t1p1_1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1); // Swap order so the oldest Polynomial will be overwritten
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(t2_2, t2p1_2); // Swap order so the oldest Polynomial will be overwritten
                phi2->next(q, dx2_gl(jj), t2_2, t2p1_2);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
}

void Wing::regularIntegralNonlinear(size_t i, size_t j, size_t n1, size_t n2, double xm, double xp, double ym, double yp)
{
    double dx = (xp - xm)/2;
    double dy = (yp - ym)/2;
    std::vector<fastgl::QuadPair> gl_x(n1), gl_y(n2);
    arma::vec x1_gl(n1, arma::fill::none), x2_gl(n2, arma::fill::none);
    for (size_t ii = 0; ii < n1; ii++)
    {
        gl_x[ii] = fastgl::GLPair(n1, ii+1);
        x1_gl(ii) = (1 - gl_x[ii].x()) * dx + xm;
    }
    for (size_t jj = 0; jj < n2; jj++)
    {
        gl_y[jj] = fastgl::GLPair(n2, jj+1);
        x2_gl(jj) = (1 - gl_y[jj].x()) * dy + ym;
    }
    arma::mat Tx_gl = Lagrange::interpolationMatrix(x1, x1_gl);
    arma::mat Ty_gl = Lagrange::interpolationMatrix(x2, x2_gl);
    arma::mat z_gl  = Lagrange::interpolation2D(Tx_gl, Ty_gl, z, x1_gl, x2_gl);
    arma::mat D1_gl = Lagrange::derivativeMatrix(x1_gl);
    arma::mat D2_gl = Lagrange::derivativeMatrix(x2_gl);
    arma::mat dz_gldx1 = D1_gl*z_gl;
    arma::mat dz_gldx2 = z_gl*D2_gl.t();
    auto [x_gl, y_gl] = Lagrange::TransfiniteQuadMap(x1_gl, x2_gl, chi);
    auto [dx_gldx1, dx_gldx2, dy_gldx1, dy_gldx2] = Lagrange::TransfiniteQuadMetrics(x1_gl, x2_gl, chi);
    arma::field<arma::mat> J_gl = {{dx_gldx1, dx_gldx2}, {dy_gldx1, dy_gldx2}};
    arma::cube e_c_gl = MetricCo(J_gl);
    arma::cube ec_gl  = MetricContra(e_c_gl);
    arma::mat e_gl = e_c_gl.slice(0)%e_c_gl.slice(2) - pow(e_c_gl.slice(1), 2);
    arma::mat sqrt_a = sqrt(e_gl%(1 + ec_gl.slice(0)%pow(dz_gldx1, 2) + 2*ec_gl.slice(1)%dz_gldx1%dz_gldx2 + ec_gl.slice(2)%pow(dz_gldx2, 2)));
    for (size_t ii = 0; ii < n1; ii++)
    {
        double  w1 = phi1->weightFunction(x1_gl(ii));
        double dw1 = phi1->weightFunctionDerivative(x1_gl(ii));
        for (size_t jj = 0; jj < n2; jj++)
        {
            arma::vec::fixed<3> r = {x_gl(ii, jj) - xC(i, j), y_gl(ii, jj) - yC(i, j), z_gl(ii, jj) - zC(i, j)};
            double r3 = pow(norm(r), 3);
            arma::vec::fixed<3> n_gl = arma::vec::fixed<3>({dy_gldx1(ii, jj)*dz_gldx2(ii, jj)-dz_gldx1(ii, jj)*dy_gldx2(ii, jj),
                                                            dz_gldx1(ii, jj)*dx_gldx2(ii, jj)-dx_gldx1(ii, jj)*dz_gldx2(ii, jj),
                                                            dx_gldx1(ii, jj)*dy_gldx2(ii, jj)-dy_gldx1(ii, jj)*dx_gldx2(ii, jj)})/sqrt_a(ii, jj);
            arma::mat::fixed<3, 2> J_red = {{dy_gldx2(ii, jj)*n_gl(2) - n_gl(1)*dz_gldx2(ii, jj),n_gl(1)*dz_gldx1(ii, jj) - dy_gldx1(ii, jj)*n_gl(2)},
                                            {n_gl(0)*dz_gldx2(ii, jj) - dx_gldx2(ii, jj)*n_gl(2),dx_gldx1(ii, jj)*n_gl(2) - n_gl(0)*dz_gldx1(ii, jj)},
                                            {dx_gldx2(ii, jj)*n_gl(1) - n_gl(0)*dy_gldx2(ii, jj),n_gl(0)*dy_gldx1(ii, jj) - dx_gldx1(ii, jj)*n_gl(1)}};
            
            double  w2 = phi2->weightFunction(x2_gl(jj));
            double dw2 = phi2->weightFunctionDerivative(x2_gl(jj));

            double t2   = phi2->constant(); // First Polynomial (j)
            double t2p1 = phi2->linear(x2_gl(jj)); // Second Polynomial (j+1 -> jp1)

            double dt2   = phi2->constantDerivative();
            double dt2p1 = phi2->linearDerivative();
            for (size_t q = 0; q < ny; q++) // Loop over Polynomial 2-direction
            {
                double  psi2 = w2*t2;
                double dpsi2 = w2*dt2 + dw2*t2;

                double t1   = phi1->constant(); // First Polynomial (i)
                double t1p1 = phi1->linear(x1_gl(ii)); // Second Polynomial (i+1 -> ip1)

                double dt1   = phi1->constantDerivative();
                double dt1p1 = phi1->linearDerivative();
                for (size_t p = 0; p < nx; p++) // Loop over Polynomial 1-direction
                {
                    double  psi1 = w1*t1;
                    double dpsi1 = w1*dt1 + dw1*t1;
                    arma::vec::fixed<2> dmudxi = {dpsi1 * psi2, psi1 * dpsi2};
                    arma::vec::fixed<3> gradmu = J_red*dmudxi;
                    arma::vec::fixed<3> gamma_gl = cross(gradmu, n_gl);
                    arma::vec q_mu = gl_x[ii].weight * gl_y[jj].weight * cross(gamma_gl, r)/r3 * dx * dy;
                    A(i+j*nx, p+q*nx) += dot(q_mu, nC.row(i+j*nx));
                    std::swap(t1, t1p1); // Swap order so the oldest Polynomial will be overwritten
                    phi1->next(p, x1_gl(ii), t1, t1p1);

                    std::swap(dt1, dt1p1);
                    phi1->nextDerivative(p, x1_gl(ii), t1, dt1, dt1p1);
                }
                std::swap(t2, t2p1); // Swap order so the oldest Polynomial will be overwritten
                phi2->next(q, x2_gl(jj), t2, t2p1);

                std::swap(dt2, dt2p1);
                phi2->nextDerivative(q, x2_gl(jj), t2, dt2, dt2p1);
            }
        }
    }
}