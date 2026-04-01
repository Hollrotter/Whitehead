#include "Structure.hpp"

void Structure::checkMesh()
{
    for (size_t i = 0; i < membranes.size(); i++)
    {
        printf("Checking mesh of membrane number %lu\n", i);
        membranes[i]->checkMesh();
    }
}

void Structure::principalStresses(const std::string nxy, const std::string n12)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->principalStresses(nxy+"_"+std::to_string(k), n12+"_"+std::to_string(k));
}

void Structure::principalStrains(const std::string vxy, const std::string gxy, const std::string g12)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->principalStrains(vxy+"_"+std::to_string(k), gxy+"_"+std::to_string(k), g12+"_"+std::to_string(k));
}

void Structure::output(const std::string filename, const Field field)
{
    for (size_t k = 0; k < membranes.size(); k++)
        membranes[k]->output(filename+"_"+std::to_string(k), field);
}