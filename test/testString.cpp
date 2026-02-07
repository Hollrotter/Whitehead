#include "String.hpp"

int main()
{
    double c = 1;
    String s(c, 101, 0.5);
    s.load(1);
    s.youngsModulus(1);
    s(Material::extensible);
    s.boundary(Location::Front, BC::Dirichlet);
    s.boundary(Location::Back,  BC::Dirichlet);
    s.linear();
    s.output("plot/Data/String/Linear");
}