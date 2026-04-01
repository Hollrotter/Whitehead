#include "VLM.hpp"

int main()
{
    switch (0)
    {
        case 0: // Rectangular Wing
        {
            double b = 10;
            double c = 1;

            arma::vec tx = arma::linspace(0, arma::datum::pi, 40);
            arma::vec ty = arma::linspace(arma::datum::pi/2, 0, 60);
            arma::vec x1 = c/2*(1 - cos(tx));
            arma::vec x2 = b/2*cos(ty);
            arma::mat x = (x1 - c/4)*arma::ones(1, x2.size());
            arma::vec y = x2;

            VLM wing(x, y);
            wing.pitch(5);
            wing(Symmetry::y);
            wing(Analysis::nonlinear);
            Camber cam([&](double xc){ return 0.1*(0.25 - pow(xc-0.5, 2)); },
                       [&](double xc){ return 0.2*(0.5-xc); });
            wing(cam);
            wing.vlm();

            double area  = 2*wing.get_area();
            arma::vec cL = wing.get_lift()   / area;
            arma::vec cM = wing.get_moment() / area/c;

            std::cout << "A  = "    << area   << '\n';
            std::cout << "cL    = " << cL.t() << '\n';
            std::cout << "cM    = " << cM.t() << '\n';

            wing.output("plot/Data/VLM/vlm");
            break;
        }
        case 1: // Rigging
        {
            double b = 15;
            double c1 = 1.5;
            double c2 = 0.9*c1;
            double alpha = 5;

            arma::vec ty = arma::linspace(0, arma::datum::pi/2, 26);
            arma::vec y  = b/2*sin(ty);
            arma::vec tx = arma::linspace(0, arma::datum::pi, 11);
            arma::mat x  = c1/2*(1-cos(tx))*(1 - (1 - c2/c1)*y/(b/2)).t();
            
            VLM wing(x, y);
            wing.pitch(alpha);
            wing.rigging({1, 0});
            wing.symmetry(Symmetry::y);
            wing.vlm();

            double     A = wing.get_area();
            arma::vec cL = wing.get_lift()   / A;
            arma::vec cM = wing.get_moment() / (A*c1);

            std::cout << "A     = " << A      << '\n';
            std::cout << "alpha = " << alpha  << '\n';
            std::cout << "cL    = " << cL.t() << '\n';
            std::cout << "cM    = " << cM.t() << '\n';

            wing.output("plot/Data/VLM/rigging");
            break;
        }
        case 2: // z0 is non-zero
        {
            arma::vec alpha = {0, 2};
            double yt = 2;
            double zt = 1;
            double cr = 1.5;
            double ct = 1.0;

            arma::mat data;
            data.load(arma::csv_name("clarky-il.csv", arma::csv_opts::no_header));
            Camber cam(Splinefit(data.col(0)/1e3, data.col(1)/1e3, 10));

            arma::vec ty = arma::linspace(0, arma::datum::pi, 21);
            arma::vec y = yt/2*(1 - cos(ty));
            arma::vec z = zt*y/yt;

            arma::vec tx = arma::linspace(0, arma::datum::pi, 21);
            arma::mat x = cr/2*(1 - cos(tx))*(1 - (1-ct/cr)*z/zt).t();

            VLM wing(x, y, z);
            wing.pitch(alpha);
            wing(Symmetry::y);
            wing(cam);
            wing.vlm();
            double     A = wing.get_area();
            arma::vec cL = wing.get_lift()   /  A;
            arma::vec cM = wing.get_moment() / (A*cr);
            auto dcLda = arma::diff(cL)/arma::diff(alpha);

            std::cout << "A     = " << A         << '\n';
            std::cout << "alpha = " << alpha.t() << '\n';
            std::cout << "cL    = " << cL.t()    << '\n';
            std::cout << "cM    = " << cM.t()    << '\n';
            std::cout << "dcLda = " << dcLda     << '\n';

            wing.output("plot/Data/VLM/z");
            break;
        }
        case 3: // Elliptical Wing
        {
            double b = 10;
            double c = 1.39;
            arma::vec alpha = {0, 5};

            arma::mat data;
            data.load(arma::csv_name("clarky-il.csv", arma::csv_opts::no_header));

            arma::vec ty = arma::linspace(arma::datum::pi/2, 0, 31);
            arma::vec x2 = 0.99*b/2*cos(ty);
            arma::vec tx = arma::linspace(0, arma::datum::pi,   21);
            arma::vec x1 = c/2*(1 - cos(tx));
            arma::mat x  = (x1 - c/4)*sqrt(1-pow(x2/(b/2), 2)).t();
            arma::vec y  = x2;

            VLM wing(x, y);
            wing.pitch(alpha);
            wing(Camber(Splinefit(data.col(0)/1e3, data.col(1)/1e3, 10)));
            wing(Symmetry::y);
            wing(Analysis::linear);
            wing.vlm();

            double A = wing.get_area();
            arma::vec cL = wing.get_lift()   /  A;
            arma::vec cM = wing.get_moment() / (A*c);
            arma::vec dcLda = 180/arma::datum::pi*arma::diff(cL)/arma::diff(alpha);

            std::cout << "A  = "      << A         << '\n';
            std::cout << "\u03B1  = " << alpha.t() << '\n';
            std::cout << "cL =  "     << cL.t()    << '\n';
            std::cout << "cM =  "     << cM.t()    << '\n';
            std::cout << "dcLda = "   << dcLda     << '\n';

            wing.output("plot/Data/VLM/EllipticalWing");
            break;
        }
    }
}