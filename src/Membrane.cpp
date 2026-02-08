#include "Membrane.hpp"

/**
 * @brief 
 * 
 * @param _Et 
 */
void Membrane::youngsModulus(const double _Et)
{
    if (_Et < 0)
    {
        std::println("Young's modulus times thickness must be positive!");
        exit(EXIT_FAILURE);
    }
    Et = _Et;
    D = Et/(1 - nu*nu);
}

/**
 * @brief 
 * 
 * @param _nu 
 */
void Membrane::poissonsRatio(const double _nu)
{
    if (_nu < 0 || _nu > 0.5)
    {
        std::println("Poisson's ratio must lie between 0 and 0.5!");
        exit(EXIT_FAILURE);
    }
    nu = _nu;
    D = Et/(1 - nu*nu);
}

void Membrane::checkMesh()
{
    printf("Checking for negative volumes...\n");
    bool negativeVolumes = false;
    arma::mat detJ = J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0);
    for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
            if (detJ(i, j) < 0)
                negativeVolumes = true;
    if (negativeVolumes == false)
        printf("No negative volumes found!\n");
    else
        printf("Negative volumes were found!\n");
}

/**
 * @brief 
 * 
 * @param f Given analytic function that defines the pressure on the nodes.
 */
void Membrane::load(const std::function<double(double, double)> f)
{
    #pragma omp parallel for
    for (size_t j = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++)
            p(i, j) = f(x(i, j), y(i, j));
}

void Membrane::inPlaneTensorial(const arma::vec _p1, const arma::vec _p2)
{
    p1 = _p1;
    p2 = _p2;
}

void Membrane::inPlaneKartesian(const arma::vec px, const arma::vec py)
{
    arma::mat det = vectorise(J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0));

    p1 = ( vectorise(J(1, 1))%px - vectorise(J(0, 1))%py)/det;
    p2 = (-vectorise(J(1, 0))%px + vectorise(J(0, 0))%py)/det;
}

void Membrane::inPlaneKartesian(const arma::mat nxx, const arma::mat nxy, const arma::mat nyy)
{
    arma::mat det2 = pow(J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0), 2);

    n11 = (pow(J(1, 1), 2)%nxx - 2*J(1, 1)%J(0, 1)%nxy + pow(J(0, 1), 2)%nyy)/det2;
    n12 =-(J(1, 1)%J(1, 0)%nxx - (J(1, 1)%J(0, 0) + J(1, 0)%J(0, 1))%nxy + J(0, 1)%J(0, 0)%nyy)/det2;
    n22 = (pow(J(1, 0), 2)%nxx - 2*J(1, 0)%J(0, 0)%nxy + pow(J(0, 0), 2)%nyy)/det2;
}

/**
 * @brief 
 * 
 * @param f Matrix containing the function values to be integrated.
 * @return double 
 */
double Membrane::integrate(arma::mat f)
{
    f.row(0).zeros();
    f.col(0).zeros();
    arma::mat Dx = D1;
    Dx.row(0).eye();
    arma::mat Dy = D2;
    Dy.row(0).eye();
    arma::mat e = e_c.slice(0)%e_c.slice(2) - e_c.slice(1)%e_c.slice(1);
    arma::mat int_f_dx   = solve(Dx, sqrt(e)%f);
    arma::mat int_f_dxdy = solve(Dy, int_f_dx.t());
    return int_f_dxdy(ny-1, nx-1) - int_f_dxdy(0, 0); // Transpose not needed because only last element needed
}

/**
 * @brief 
 * 
 * @param field Defines the field that will be integrated.
 * @return double 
 */
double Membrane::integrate(const Field field)
{
    std::unique_ptr<TensorField> f = setField(field);
    arma::mat F = *f;
    f.release();
    return integrate(F);
}

/**
 * @brief 
 * 
 * @return double 
 */
double Membrane::elasticPotential()
{
    return integrate(gamma_11%n11 + 2*gamma_12%n12 + gamma_22%n22)/2;
}

std::pair<arma::mat, arma::mat> Membrane::kartesianDisplacements()
{
    arma::mat det = J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0);

    arma::mat vx = ( J(1, 1)%v1 - J(1, 0)%v2)/det;
    arma::mat vy = (-J(0, 1)%v1 + J(0, 0)%v2)/det;

    return std::tie(vx, vy);
}

void Membrane::principalStresses(const std::string filenameKartesianStresses, const std::string filenamePrincipalStresses)
{
    arma::mat sigma_x = pow(J(0, 0), 2)%n11 +                   2*J(0, 0)%J(0, 1)%n12 + pow(J(0, 1), 2)%n22;
    arma::mat tau_xy  = J(0, 0)%J(1, 0)%n11 + (J(0, 0)%J(1, 1) + J(0, 1)%J(1, 0))%n12 + J(0, 1)%J(1, 1)%n22;
    arma::mat sigma_y = pow(J(1, 0), 2)%n11 +                   2*J(1, 0)%J(1, 1)%n12 + pow(J(1, 1), 2)%n22;

    arma::mat sigma_1 = (sigma_x + sigma_y)/2 + sqrt((sigma_x - sigma_y)%(sigma_x - sigma_y)/4 + tau_xy%tau_xy);
    arma::mat sigma_2 = (sigma_x + sigma_y)/2 - sqrt((sigma_x - sigma_y)%(sigma_x - sigma_y)/4 + tau_xy%tau_xy);
    arma::mat tau_12  = (sigma_1 - sigma_2)/2;
    std::ofstream file_k(filenameKartesianStresses);
    std::ofstream file_p(filenamePrincipalStresses);
    for (size_t i = 0; i < nx; i++, file_k<<'\n', file_p<<'\n')
        for (size_t j = 0; j < ny; j++, file_k<<'\n', file_p<<'\n')
        {
            file_k << x(i, j) << ' ' << y(i, j) << ' ' << sigma_x(i, j) << ' ' << sigma_y(i, j) << ' ' << tau_xy(i, j);
            file_p << x(i, j) << ' ' << y(i, j) << ' ' << sigma_1(i, j) << ' ' << sigma_2(i, j) << ' ' << tau_12(i, j);
        }
    file_k.close();
    file_p.close();
}

void Membrane::principalStrains(const std::string filenameKartesianDeformations, const std::string filenameKartesianStrains, const std::string filenamePrincipalStrains)
{
    arma::mat det2 = pow(J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0), 2);

    auto [vx, vy] = kartesianDisplacements();

    arma::mat epsilon_x = (pow(J(1, 1), 2)%gamma_11 - 2*J(1, 1)%J(1, 0)%gamma_12  + pow(J(1, 0), 2)%gamma_22)/det2;
    arma::mat gamma_xy  =-(J(1, 1)%J(0, 1)%gamma_11 - (J(1, 1)%J(0, 0) + J(1, 0)%J(0, 1))%gamma_12 + J(1, 0)%J(0, 0)%gamma_22)/det2;
    arma::mat epsilon_y = (pow(J(0, 1), 2)%gamma_11 - 2*J(0, 1)%J(0, 0)%gamma_12  + pow(J(0, 0), 2)%gamma_22)/det2;
    
    arma::mat epsilon_1 = (epsilon_x + epsilon_y)/2 + sqrt((epsilon_x - epsilon_y)%(epsilon_x - epsilon_y)/4 + gamma_xy%gamma_xy);
    arma::mat epsilon_2 = (epsilon_x + epsilon_y)/2 - sqrt((epsilon_x - epsilon_y)%(epsilon_x - epsilon_y)/4 + gamma_xy%gamma_xy);
    arma::mat gamma_12  = (epsilon_1 - epsilon_2)/2;
    std::ofstream file_v(filenameKartesianDeformations);
    std::ofstream file_k(filenameKartesianStrains);
    std::ofstream file_p(filenamePrincipalStrains);
    for (size_t i = 0; i < nx; i++, file_v<<'\n', file_k<<'\n', file_p<<'\n')
        for (size_t j = 0; j < ny; j++, file_v<<'\n', file_k<<'\n', file_p<<'\n')
        {
            file_v << x(i, j) << ' ' << y(i, j) << ' ' << vx(i, j) << ' ' << vy(i, j);
            file_k << x(i, j) << ' ' << y(i, j) << ' ' << epsilon_x(i, j) << ' ' << epsilon_y(i, j) << ' ' << gamma_xy(i, j);
            file_p << x(i, j) << ' ' << y(i, j) << ' ' << epsilon_1(i, j) << ' ' << epsilon_2(i, j) << ' ' << gamma_12(i, j);
        }
    file_v.close();
    file_k.close();
    file_p.close();
}

/**
 * @brief 
 * 
 * @param filename Path of the file where the data will be output to.
 * @param field Field to be output.
 */
void Membrane::output(const std::string filename, const Field field)
{
    std::unique_ptr<TensorField> f = setField(field);
    std::ofstream file(filename);
    for (size_t i = 0; i < nx; i++, file<<'\n')
        for (size_t j = 0; j < ny; j++, file<<'\n')
            file << x(i, j) << ' ' << y(i, j) << ' ' << f->operator()(i, j);
    file.close();
    f.release();
}

double Membrane::operator()(const size_t i, const size_t j, const Field field)
{
    std::unique_ptr<TensorField> f = setField(field);
    double F = f->operator()(i, j);
    f.release();
    return F;
}

TensorField Membrane::operator()(const Field field)
{
    std::unique_ptr<TensorField> f = setField(field);
    TensorField F = *f;
    f.release();
    return F;
}

/**
 * @brief 
 * 
 * @param field Field the pointer will point at.
 * @return std::unique_ptr<TensorField> 
 */
std::unique_ptr<TensorField> Membrane::setField(const Field field)
{
    std::unique_ptr<TensorField> f;
    switch(field)
    {
        case Field::z:
            f.reset(&z);
            return f;
        case Field::v1:
            f.reset(&v1);
            return f;
        case Field::v2:
            f.reset(&v2);
            return f;
        case Field::n11:
            f.reset(&n11);
            return f;
        case Field::n12:
            f.reset(&n12);
            return f;
        case Field::n22:
            f.reset(&n22);
            return f;
        default:
            std::println("Not a valid field for the Membrane structure!");
            f.release();
            exit(EXIT_FAILURE);
    }
    std::unreachable();
}

void Membrane::operator()(const arma::mat xNew, const arma::mat yNew)
{
    size_t nxNew = xNew.n_rows;
    size_t nyNew = yNew.n_cols;
    nxy = nxNew*nyNew;

    arma::vec x1New = Chebyshev::gaussLobatto(nxNew);
    arma::vec x2New = Chebyshev::gaussLobatto(nyNew);

    arma::mat T1 = Lagrange::interpolationMatrix(x1, x1New);
    arma::mat T2 = Lagrange::interpolationMatrix(x2, x2New);

    z = Lagrange::interpolation2D(T1, T2, z, x1New, x2New);

    arma::mat det = J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0);

    arma::mat vx = ( J(1, 1)%v1 - J(1, 0)%v2)/det;
    arma::mat vy = (-J(0, 1)%v1 + J(0, 0)%v2)/det;

    vx  = Lagrange::interpolation2D(T1, T2,  vx, x1New, x2New);
    vy  = Lagrange::interpolation2D(T1, T2,  vy, x1New, x2New);

    p   = Lagrange::interpolation2D(T1, T2, p, x1New, x2New);
    p1  = arma::vectorise(Lagrange::interpolation2D(T1, T2, arma::reshape(p1, nx, ny), x1New, x2New));
    p2  = arma::vectorise(Lagrange::interpolation2D(T1, T2, arma::reshape(p2, nx, ny), x1New, x2New));

    nx = nxNew;
    ny = nyNew;
    x1 = x1New;
    x2 = x2New;
    x  =  xNew;
    y  =  yNew;

    D1  = Chebyshev::derivativeMatrix(x1, Derivative::first);
    D2  = Chebyshev::derivativeMatrix(x2, Derivative::first);
    D11 = Chebyshev::derivativeMatrix(x1, Derivative::second);
    D22 = Chebyshev::derivativeMatrix(x2, Derivative::second);
    DD1   = ddx();
    DD2   = ddy();
    D2D11 = d2dx2();
    D2D12 = d2dxdy();
    D2D22 = d2dy2();
    J = Jacobian(x, y, D1, D2);
    e_c = MetricCo(J);
    ec  = MetricContra(e_c);
    e = e_11()%e_22() - pow(e_12(), 2);
    gam = Christoffel(e_c, ec, D1, D2);
    Psi1111 = 2*(g111()%g111() + g211()%g112()) - D1*g111();
    Psi1112 = 2*(g111()%g211() + g211()%g212()) - D1*g211();
    Psi2221 = 2*(g122()%g112() + g222()%g122()) - g122()*D2.t();
    Psi2222 = 2*(g122()%g212() + g222()%g222()) - g222()*D2.t();
    Psi1121 = 2*(g111()%g112() + g112()%g212()) - g111()*D2.t();
    Psi1122 = 2*(g112()%g211() + g212()%g212()) - g211()*D2.t();
    Psi2121 = 2*(g112()%g112() + g212()%g122()) - D1*g122();
    Psi2122 = 2*(g112()%g212() + g212()%g222()) - D1*g222();
    Psi1221 = g122()%g111() + g222()%g112() + g112()%g112() + g212()%g122() - g112()*D2.t();
    Psi1222 = g122()%g211() + g112()%g212() + 2*g212()%g222() - g212()*D2.t();
    Psi2111 = g111()%g212() + g211()%g222() + g112()%g211() + g212()%g212() - D1*g212();
    Psi2112 = 2*g111()%g112() + g211()%g122() + g212()%g112() - D1*g112();
    H1111 = e11()%e11();
    H1112 = e11()%e12();
    H1222 = e12()%e22();
    H1122 = e12()%e12() + nu/e;
    H1221 = (e11()%e22() + e12()%e12() - nu/e)/2;
    H2222 = e22()%e22();

    h_1s1_east  = sqrt(e11().row(nx-1));
    h_11s_east  = 1/h_1s1_east;
    h_22s_east  = sqrt(e_22().row(nx-1));
    h_2s2_east  = 1/h_22s_east;
    h_12s_east  = e_12().row(nx-1)/sqrt(e_22().row(nx-1));
    h_1s2_east  = h_11s_east%e12().row(nx-1);
    h_1s1_west  = sqrt(e11().row(0));
    h_11s_west  = 1/h_1s1_west;
    h_22s_west  = sqrt(e_22().row(0));
    h_2s2_west  = 1/h_22s_west;
    h_12s_west  = e_12().row(0)/sqrt(e_22().row(0));
    h_1s2_west  = h_11s_west%e12().row(0);
    h_11s_south = sqrt(e_11().col(0));
    h_1s1_south = 1/h_11s_south;
    h_2s2_south = sqrt(e22().col(0));
    h_22s_south = 1/h_2s2_south;
    h_21s_south = e_12().col(0)/sqrt(e_11().col(0));
    h_2s1_south = h_22s_south%e12().col(0);
    h_11s_north = sqrt(e_11().col(ny-1));
    h_1s1_north = 1/h_11s_north;
    h_2s2_north = sqrt(e22().col(ny-1));
    h_22s_north = 1/h_2s2_north;
    h_21s_north = e_12().col(ny-1)/sqrt(e_11().col(ny-1));
    h_2s1_north = h_22s_north%e12().col(ny-1);

    v1 = J(0, 0)%vx + J(1, 0)%vy;
    v2 = J(0, 1)%vx + J(1, 1)%vy;

    v1__1 = D1*v1     - g111()%v1 - g211()%v2;
    v1__2 = v1*D2.t() - g112()%v1 - g212()%v2;
    v2__1 = D1*v2     - g112()%v1 - g212()%v2;
    v2__2 = v2*D2.t() - g122()%v1 - g222()%v2;

    if (analysis == Analysis::nonlinear)
    {
        arma::mat z_1  = D1*z;
        arma::mat z_2  = z*D2.t();
        gamma_11 =  v1__1         + (z_1%z_1 + e11()%v1__1%v1__1 + 2*e12()%v1__1%v2__1 + e22()%v2__1%v2__1)/2;
        gamma_12 = (v1__2 + v2__1 +  z_1%z_2 + e11()%v1__1%v1__2 + e12()%(v1__1%v2__2 + v2__1%v1__2) + e22()%v2__1%v2__2)/2;
        gamma_22 =  v2__2         + (z_2%z_2 + e11()%v1__2%v1__2 + 2*e12()%v1__2%v2__2 + e22()%v2__2%v2__2)/2;
    }
    else
    {
        gamma_11 = v1__1;
        gamma_12 = (v1__2 + v2__1)/2;
        gamma_22 = v2__2;
    }

    n11 = D*(H1111%gamma_11 + 2*H1112%gamma_12 + H1122%gamma_22);
    n12 = D*(H1112%gamma_11 + 2*H1221%gamma_12 + H1222%gamma_22);
    n22 = D*(H1122%gamma_11 + 2*H1222%gamma_12 + H2222%gamma_22);

    S = arma::zeros(nxy, nxy);
    L = arma::zeros(nxy, nxy);
    U = arma::zeros(nxy, nxy);
    P = arma::zeros(nxy, nxy);
    b = arma::zeros(nxy);

    z.initBC(nxNew, nyNew);
    v1.initBC(nxNew, nyNew);
    v2.initBC(nxNew, nyNew);
    n11.initBC(nxNew, nyNew);
    n12.initBC(nxNew, nyNew);
    n22.initBC(nxNew, nyNew);
}

void Membrane::operator()(const std::array<Lagrange::CurveInterpolant*, 4> _chi)
{
    chi = _chi;
    size_t nxNew = chi[1]->getNodes().size();
    size_t nyNew = chi[0]->getNodes().size();
    nxy = nxNew*nyNew;

    arma::vec x1New = Chebyshev::gaussLobatto(nxNew);
    arma::vec x2New = Chebyshev::gaussLobatto(nyNew);

    arma::mat T1 = Lagrange::interpolationMatrix(x1, x1New);
    arma::mat T2 = Lagrange::interpolationMatrix(x2, x2New);

    auto [xNew, yNew] = Lagrange::TransfiniteQuadMap(chi);

    z   = Lagrange::interpolation2D(T1, T2,   z, x1New, x2New);

    arma::mat det = J(0, 0)%J(1, 1) - J(0, 1)%J(1, 0);

    arma::mat vx = ( J(1, 1)%v1 - J(1, 0)%v2)/det;
    arma::mat vy = (-J(0, 1)%v1 + J(0, 0)%v2)/det;

    vx  = Lagrange::interpolation2D(T1, T2,  vx, x1New, x2New);
    vy  = Lagrange::interpolation2D(T1, T2,  vy, x1New, x2New);

    p   = Lagrange::interpolation2D(T1, T2, p, x1New, x2New);
    p1  = arma::vectorise(Lagrange::interpolation2D(T1, T2, arma::reshape(p1, nx, ny), x1New, x2New));
    p2  = arma::vectorise(Lagrange::interpolation2D(T1, T2, arma::reshape(p2, nx, ny), x1New, x2New));

    nx = nxNew;
    ny = nyNew;
    x1 = x1New;
    x2 = x2New;
    x  =  xNew;
    y  =  yNew;

    D1  = Chebyshev::derivativeMatrix(x1, Derivative::first);
    D2  = Chebyshev::derivativeMatrix(x2, Derivative::first);
    D11 = Chebyshev::derivativeMatrix(x1, Derivative::second);
    D22 = Chebyshev::derivativeMatrix(x2, Derivative::second);
    DD1   = ddx();
    DD2   = ddy();
    D2D11 = d2dx2();
    D2D12 = d2dxdy();
    D2D22 = d2dy2();
    J = Jacobian(xNew, yNew, D1, D2);
    e_c = MetricCo(J);
    ec  = MetricContra(e_c);
    e = e_11()%e_22() - pow(e_12(), 2);
    gam = Christoffel(e_c, ec, D1, D2);
    Psi1111 = 2*(g111()%g111() + g211()%g112()) - D1*g111();
    Psi1112 = 2*(g111()%g211() + g211()%g212()) - D1*g211();
    Psi2221 = 2*(g122()%g112() + g222()%g122()) - g122()*D2.t();
    Psi2222 = 2*(g122()%g212() + g222()%g222()) - g222()*D2.t();
    Psi1121 = 2*(g111()%g112() + g112()%g212()) - g111()*D2.t();
    Psi1122 = 2*(g112()%g211() + g212()%g212()) - g211()*D2.t();
    Psi2121 = 2*(g112()%g112() + g212()%g122()) - D1*g122();
    Psi2122 = 2*(g112()%g212() + g212()%g222()) - D1*g222();
    Psi1221 = g122()%g111() + g222()%g112() + g112()%g112() + g212()%g122() - g112()*D2.t();
    Psi1222 = g122()%g211() + g112()%g212() + 2*g212()%g222() - g212()*D2.t();
    Psi2111 = g111()%g212() + g211()%g222() + g112()%g211() + g212()%g212() - D1*g212();
    Psi2112 = 2*g111()%g112() + g211()%g122() + g212()%g112() - D1*g112();
    H1111 = e11()%e11();
    H1112 = e11()%e12();
    H1222 = e12()%e22();
    H1122 = e12()%e12() + nu/e;
    H1221 = (e11()%e22() + e12()%e12() - nu/e)/2;
    H2222 = e22()%e22();

    h_1s1_east = sqrt(e11().row(nx-1));
    h_11s_east = 1/h_1s1_east;
    h_22s_east = sqrt(e_22().row(nx-1));
    h_2s2_east = 1/h_22s_east;
    h_12s_east = e_12().row(nx-1)/sqrt(e_22().row(nx-1));
    h_1s2_east = h_11s_east%e12().row(nx-1);
    h_1s1_west = sqrt(e11().row(0));
    h_11s_west = 1/h_1s1_west;
    h_22s_west = sqrt(e_22().row(0));
    h_2s2_west = 1/h_22s_west;
    h_12s_west = e_12().row(0)/sqrt(e_22().row(0));
    h_1s2_west = h_11s_west%e12().row(0);
    h_11s_south = sqrt(e_11().col(0));
    h_1s1_south = 1/h_11s_south;
    h_2s2_south = sqrt(e22().col(0));
    h_22s_south = 1/h_2s2_south;
    h_21s_south = e_12().col(0)/sqrt(e_11().col(0));
    h_2s1_south = h_22s_south%e12().col(0);
    h_11s_north = sqrt(e_11().col(ny-1));
    h_1s1_north = 1/h_11s_north;
    h_2s2_north = sqrt(e22().col(ny-1));
    h_22s_north = 1/h_2s2_north;
    h_21s_north = e_12().col(ny-1)/sqrt(e_11().col(ny-1));
    h_2s1_north = h_22s_north%e12().col(ny-1);
    
    v1 = J(0, 0)%reshape(vx, nx, ny) + J(1, 0)%reshape(vy, nx, ny);
    v2 = J(0, 1)%reshape(vx, nx, ny) + J(1, 1)%reshape(vy, nx, ny);

    v1__1 = D1*v1     - g111()%v1 - g211()%v2;
    v1__2 = v1*D2.t() - g112()%v1 - g212()%v2;
    v2__1 = D1*v2     - g112()%v1 - g212()%v2;
    v2__2 = v2*D2.t() - g122()%v1 - g222()%v2;

    if (analysis == Analysis::nonlinear)
    {
        arma::mat z_1  = D1*z;
        arma::mat z_2  = z*D2.t();
        gamma_11 =  v1__1         + (z_1%z_1 + e11()%v1__1%v1__1 + 2*e12()%v1__1%v2__1 + e22()%v2__1%v2__1)/2;
        gamma_12 = (v1__2 + v2__1 +  z_1%z_2 + e11()%v1__1%v1__2 + e12()%(v1__1%v2__2 + v2__1%v1__2) + e22()%v2__1%v2__2)/2;
        gamma_22 =  v2__2         + (z_2%z_2 + e11()%v1__2%v1__2 + 2*e12()%v1__2%v2__2 + e22()%v2__2%v2__2)/2;
    }
    else
    {
        gamma_11 = v1__1;
        gamma_12 = (v1__2 + v2__1)/2;
        gamma_22 = v2__2;
    }

    n11 = D*(H1111%gamma_11 + 2*H1112%gamma_12 + H1122%gamma_22);
    n12 = D*(H1112%gamma_11 + 2*H1221%gamma_12 + H1222%gamma_22);
    n22 = D*(H1122%gamma_11 + 2*H1222%gamma_12 + H2222%gamma_22);

    S = arma::zeros(nxy, nxy);
    L = arma::zeros(nxy, nxy);
    U = arma::zeros(nxy, nxy);
    P = arma::zeros(nxy, nxy);
    b = arma::zeros(nxy);

    z.initBC(nxNew, nyNew);
    v1.initBC(nxNew, nyNew);
    v2.initBC(nxNew, nyNew);
    n11.initBC(nxNew, nyNew);
    n12.initBC(nxNew, nyNew);
    n22.initBC(nxNew, nyNew);
}