#include "Structure.hpp"

int main()
{
    size_t n1 = 11;
    size_t n2 = 11;
    size_t n3 = 11;

    double a = 1;
    double Et = 1000;
    double nu = 0.3;
    double D = Et/(1 - nu*nu);

    Point p1(0, 0);
    Point p2(a/2, 0);
    Point p3(a/(6-sqrt(12)), a/(6-sqrt(12)));
    Point p4(0, a/2);
    Point p5(a, 0);
    Point p6(0, a);
    Point p7(a/sqrt(2), a/sqrt(2));

    Lagrange::CurveInterpolant chi1(p1, p2, n1);
    Lagrange::CurveInterpolant chi2(p2, p3, n2);
    Lagrange::CurveInterpolant chi3(p4, p3, n1);
    Lagrange::CurveInterpolant chi4(p1, p4, n2);
    Lagrange::CurveInterpolant chi5(p2, p5, n3);
    Lagrange::CurveInterpolant chi6(p4, p6, n3);
    Lagrange::CurveInterpolant chi7(p3, p7, n3);
    Lagrange::CurveInterpolant chi8(p5, p7, p1, a, n2);
    Lagrange::CurveInterpolant chi9(p6, p7, p1, a, n1);

    Membrane m1({&chi1, &chi2, &chi3, &chi4});
    Membrane m2({&chi5, &chi8, &chi7, &chi2});
    Membrane m3({&chi3, &chi7, &chi9, &chi6});

    m1.boundary(Field::z,   &chi1, BC::Neumann);
    m1.boundary(Field::z,   &chi4, BC::Neumann);
    m2.boundary(Field::z,   &chi5, BC::Neumann);
    m2.boundary(Field::z,   &chi8, BC::Dirichlet);
    m3.boundary(Field::z,   &chi6, BC::Neumann);
    m3.boundary(Field::z,   &chi9, BC::Dirichlet);

    m1.boundary(Field::v1,  &chi1, BC::Neumann);
    m1.boundary(Field::v1,  &chi4, BC::Dirichlet);
    m2.boundary(Field::v1,  &chi5, BC::Neumann);
    m3.boundary(Field::v1,  &chi6, BC::Dirichlet);

    m1.boundary(Field::v2,  &chi1, BC::Dirichlet);
    m1.boundary(Field::v2,  &chi4, BC::Neumann);
    m2.boundary(Field::v2,  &chi5, BC::Dirichlet);
    m3.boundary(Field::v2,  &chi6, BC::Neumann);

    m2.boundary(Field::n11, &chi8, BC::Dirichlet, 0.005*D);
    m2.boundary(Field::n12, &chi8, BC::Dirichlet);
    m3.boundary(Field::n12, &chi9, BC::Dirichlet);
    m3.boundary(Field::n22, &chi9, BC::Dirichlet, 0.005*D);

    Structure s({&m1, &m2, &m3});

    arma::vec x = join_vert(m1.X().col(0), m2.X().col(0));
    arma::vec z = join_vert(m1.Z().col(0), m2.Z().col(0));

    std::ofstream file0("plot/Data/CircularMembrane_0.txt");
    for (size_t i = 0; i < z.size(); i++)
        file0 << x(i) << ' ' << z(i) << '\n';
    file0.close();

    s.output("plot/Data/CircularMembrane_z", Field::z);

    arma::mat J11_1 = m1.J11();
    arma::mat J12_1 = m1.J12();
    arma::mat J21_1 = m1.J21();
    arma::mat J22_1 = m1.J22();
    arma::mat J11_2 = m2.J11();
    arma::mat J12_2 = m2.J12();
    arma::mat J21_2 = m2.J21();
    arma::mat J22_2 = m2.J22();

    s.youngsModulus(Et);
    s.poissonsRatio(nu);
    s.setIterations(100);
    s.planeStrain();
    s.load(0.01*D/a);
    s.semilinear();

    m2.boundary(Field::n11, &chi8, BC::None);
    m2.boundary(Field::n12, &chi8, BC::None);
    m3.boundary(Field::n12, &chi9, BC::None);
    m3.boundary(Field::n22, &chi9, BC::None);

    m2.boundary(Field::v1, &chi8, BC::Dirichlet);
    m2.boundary(Field::v2, &chi8, BC::Dirichlet);
    m3.boundary(Field::v1, &chi9, BC::Dirichlet);
    m3.boundary(Field::v2, &chi9, BC::Dirichlet);

    s.setIterations(100);

    std::cout << "kappa = 0.01\n";
    s.load(0.01*D/a);
    s.nonlinear();

    arma::vec kappa = {0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};

    for (size_t k = 0; k < kappa.size(); k++)
    {
        std::cout << "kappa = " << kappa(k) << '\n';
        s.load(kappa(k)*D/a);
        s.nonlinear();

        arma::mat vx_1 = (J22_1%m1.V1() - J21_1%m1.V2())/(J11_1%J22_1 - J12_1%J21_1);
        arma::mat vx_2 = (J22_2%m2.V1() - J21_2%m2.V2())/(J11_2%J22_2 - J12_2%J21_2);

        arma::vec x = join_vert(m1.X().col(0)+vx_1.col(0), m2.X().col(0)+vx_2.col(0));
        arma::vec z = join_vert(m1.Z().col(0), m2.Z().col(0));

        std::ofstream file("plot/Data/CircularMembrane_"+std::to_string(k+1)+".txt");
        for (size_t i = 0; i < z.size(); i++)
            file << x(i) << ' ' << z(i) << '\n';
        file.close();
    }
}