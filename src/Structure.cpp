#include "Structure.hpp"

void Structure::checkMesh()
{
    for (size_t i = 0; i < membranes.size(); i++)
    {
        printf("Checking mesh of membrane number %lu\n", i);
        membranes[i]->checkMesh();
    }
}

void Structure::principalStresses(const std::string filenameKartesianStresses, const std::string filenamePrincipalStresses)
{
    for (size_t k = 0; k < membranes.size(); k++)
    {
        std::string file_k1 = filenameKartesianStresses;
        std::string file_k2 = filenamePrincipalStresses;
        file_k1.append("_");
        file_k2.append("_");
        file_k1.append(std::to_string(k));
        file_k2.append(std::to_string(k));
        membranes[k]->principalStresses(file_k1, file_k2);
    }
}

void Structure::principalStrains(const std::string filenameKartesianDeformations, const std::string filenameKartesianStrains, const std::string filenamePrincipalStrains)
{
    for (size_t k = 0; k < membranes.size(); k++)
    {
        std::string file_k1 = filenameKartesianDeformations;
        std::string file_k2 = filenameKartesianStrains;
        std::string file_k3 = filenamePrincipalStrains;
        file_k1.append("_");
        file_k2.append("_");
        file_k3.append("_");
        file_k1.append(std::to_string(k));
        file_k2.append(std::to_string(k));
        file_k3.append(std::to_string(k));
        membranes[k]->principalStrains(file_k1, file_k2, file_k3);
    }
}

void Structure::output(const std::string filename, const Field field)
{
    for (size_t k = 0; k < membranes.size(); k++)
    {
        std::string file_k = filename;
        file_k.append("_");
        file_k.append(std::to_string(k));
        membranes[k]->output(file_k, field);
    }
}