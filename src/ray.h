#ifndef RAY_H
#define RAY_H

#include <LinearMath/btVector3.h>

#include <units/units.h>
#include "mesh.h"
#include <iostream>
#include <fstream>
// Objective: define the characteristcs of the ray should have, like what is ray, segment and hit boundary
// the ray function is the basis of combining with physics bullet to implement ray tracing algorithm
// Note: try to interpert the each parameter of struct ray and struct segment, which are important to the following applcation
 
// class material;
class mesh;

namespace ray_physics {

struct ray
{
    btVector3 from, direction;
    size_t depth;
    const material & media;
    const material * media_outside = nullptr; // keep track of media outside of vascularities, so we can switch back to it when the ray gets out of them
    float intensity, frequency;
    units::length::millimeter_t distance_traveled; // [mm]
    unsigned short parent_collision; // position in collision vector

    static constexpr size_t max_depth = 12;
    static constexpr float intensity_epsilon = 0.1; //TODO: test
    bool null = false;
};

struct segment
{
    btVector3 from, to, direction;
    float reflected_intensity; // reflected back to the transducer, at the end of the segment
    float initial_intensity, attenuation;

    units::length::millimeter_t distance_traveled; // traveled from the transducer to the beginning of the segment
    const material & media;
};

struct collision
{
    btVector3 position;
    unsigned short parent_collision; // position in collision vector
};

struct hit_result { float reflected_intensity; ray reflection, refraction; };

hit_result hit_boundary(ray & r, const btVector3 &hit_point, const btVector3 & surface_normal, const mesh & collided_mesh);

// Advance through homogeneous media and decrease intensity accordingly
units::length::millimeter_t segment_length_in_mm(const btVector3 & v1, const btVector3 & v2);
void travel(ray & r, units::length::millimeter_t mm);

bool should_travel(const ray & r);

float max_ray_length(const ray & r);

btVector3 snells_law(const btVector3 & ray_direction, const btVector3 & surface_normal, float incidence_angle, float refraction_angle, float refr_ratio);

/**
 * Intensity of the reflected ray.
 * IMPORTANT: This it NOT the intensity reflected back to the transducer.
 */
float reflection_intensity(const float intensity_in, const float media_1, const float incidence_angle, const float media_2, const float refracted_angle);

// Intensity reflected back to the transducer when a ray passes through an interface.
float reflected_intensity(const float ray_intensity, const float incidence_angle, const material & ray_media, const material & colliding_media);

btVector3 random_unit_vector(btVector3 v, float cos_theta);

float power_cosine_variate(int v);

void writeTextToFile(float index, float text, const char* filename);

} // end ray_physics

#endif // RAY_H
