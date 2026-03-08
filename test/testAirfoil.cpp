#include "Airfoil.hpp"

int main()
{
    switch (0)
    {
        case 0: // Flat plate
        {
            double c = 2;
            size_t N = 10;

            Airfoil airfoil(c, N);
            airfoil.pitch(5);
            airfoil.linear();
            airfoil.output("plot/Data/Airfoil/flat");

            std::cout << "cL = " << airfoil.get_lift().t()   << std::endl;
            std::cout << "cM = " << airfoil.get_moment().t() << std::endl;
            break;
        }
    }
}