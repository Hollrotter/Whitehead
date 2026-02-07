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