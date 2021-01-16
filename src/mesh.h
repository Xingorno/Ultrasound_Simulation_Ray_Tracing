#ifndef MESH_H
#define MESH_H

#include <string>
#include <array>

// Objective: define the mesh model and what parameters/attributes should be included
// Note: how to assign the parameters/attributes to the different parts of objects
struct material
{
    float impedance, attenuation, mu0, mu1, sigma, specularity,/*specular(+infinite) -> diffuse (0)*/ roughness;
};

struct mesh
{
    const std::string filename;
    const bool is_rigid, is_vascular;
    const std::array<float,3> deltas;
    const bool outside_normals;
    const material & material_inside;
    const material & material_outside;
};

#endif // MESH_H
