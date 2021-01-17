#include "ray.h"

#include <cmath>
#include <iostream>
#include <random>
#include "mesh.h"

using namespace ray_physics;
using namespace std;
float index_num= 1.0;
ray_physics::hit_result ray_physics::hit_boundary(const ray & r, const btVector3 & hit_point, const btVector3 & surface_normal, const mesh & collided_mesh)
{
    // btScalar incidence_angle = r.direction.dot(-surface_normal); // cos theta_1

    // TODO: this logic can probably be simpler
    const material * material_after_vascularities = nullptr;
    // const auto & material_after_collision = [&r, &collided_mesh, &material_after_vascularities]() -> const material &
    // {
    //     if (r.media_outside) // if we are in a vessel
    //     {
    //         if (collided_mesh.is_vascular) // and we collided against a vessel (assuming same vessel, so we're getting out of it)
    //         {
    //             material_after_vascularities = nullptr;
    //             return *r.media_outside; // we are going back to the stored media
    //         }
    //         // UNDONE
    //         else // we are still inside the vessel but went out of the surrounding organ
    //         {
    //             // update the surrounding tissue
    //             material_after_vascularities = r.media_outside == &collided_mesh.material_inside ? &collided_mesh.material_outside : &collided_mesh.material_inside;
    //             // but we remain in the same media, i.e. the vessel
    //             return r.media;
    //         }
    //     }
    //     else // we are not in a vessel
    //     {
    //         if (collided_mesh.is_vascular) // and we collided with a vessel
    //         {
    //             // update the surrounding tissue
    //             material_after_vascularities = &r.media; // we will come back to this tissue after getting out of the vessel
    //             return collided_mesh.material_inside;
    //         }
    //         else // and we collided with a regular organ
    //         {
    //             material_after_vascularities = nullptr;
    //             return &r.media == &collided_mesh.material_inside ? collided_mesh.material_outside : collided_mesh.material_inside;
    //         }
    //     }
    // }();
    // Given this 
    const auto & material_after_collision = collided_mesh.material_outside;
    
    //power-cosine distribution random normal
    float random_angle = power_cosine_variate(material_after_collision.roughness); // par: s shininess : ( 0 = diffuse; inf = specular)
    btVector3 random_normal = random_unit_vector(surface_normal, random_angle);
    random_normal.safeNormalize();
    //Replace the surface normal with the new random normal
    btScalar incidence_angle = r.direction.dot(-random_normal); // cos theta_1
    // Question: Why?
    if (incidence_angle < 0)
    {
        // incidence_angle = r.direction.dot(surface_normal);
        incidence_angle = r.direction.dot(random_normal);
        random_normal = -random_normal;
    }

    const float refr_ratio = r.media.impedance / material_after_collision.impedance;
    
    // // Refraction angle
    // float refraction_angle = 1 - refr_ratio*refr_ratio * (1 - incidence_angle*incidence_angle);
    
    // refraction_angle = std::sqrt(refraction_angle);

    // // Note: modify the refraction angle due to the strong refraction effect
    // float theta_incidence_rad = acos(incidence_angle);
    // float theta_refraction_rad = acos(refraction_angle);
    // refraction_angle = cos((theta_refraction_rad-theta_incidence_rad)*0/5 + theta_incidence_rad); 
    // const bool total_internal_reflection = refraction_angle < 0;
    // // Refraction direction
    // const auto refraction_direction = snells_law(r.direction, random_normal, incidence_angle, refraction_angle, refr_ratio);
    
    /* Orthognal incidence */
    float refraction_angle = incidence_angle;
    const bool total_internal_reflection = refraction_angle < 0;
    const auto refraction_direction = r.direction;
    
    // Reflection direction
    
    btVector3 reflection_direction1 = r.direction + 2*incidence_angle * random_normal;
    const btVector3 reflection_direction = reflection_direction1.safeNormalize();


    // Reflection intensity
    const auto intensity_refl = total_internal_reflection ?
                                    r.intensity :
                                    reflection_intensity(r.intensity,
                                        r.media.impedance, incidence_angle,
                                        material_after_collision.impedance, refraction_angle);
    
    // Refraction intensity
    const auto intensity_refr = r.intensity - intensity_refl;

    // Eq. 10 in Burger13
    const float back_to_transducer_intensity = reflected_intensity(r.intensity, incidence_angle, r.media, material_after_collision);

    // Add two more rays to the stack

    // units::length::millimeter_t temp = segment_length_in_mm(r.from, hit_point);
    // units::length::millimeter_t distanceTraved_new = temp + r.distance_traveled;

    // ray refraction_ray { hit_point, refraction_direction, r.depth+1, material_after_collision, material_after_vascularities, intensity_refr > ray::intensity_epsilon ? intensity_refr : 0.0f, r.frequency, distanceTraved_new, 0 };

    // ray reflection_ray { hit_point, reflection_direction, r.depth+1, r.media, r.media_outside, intensity_refl > ray::intensity_epsilon ? intensity_refl : 0.0f, r.frequency, distanceTraved_new, 0 };

    ray refraction_ray { hit_point, refraction_direction, r.depth+1, material_after_collision, material_after_vascularities, intensity_refr > ray::intensity_epsilon ? intensity_refr : 0.0f, r.frequency, r.distance_traveled, 0 };

    ray reflection_ray { hit_point, reflection_direction, r.depth+1, r.media, r.media_outside, intensity_refl > ray::intensity_epsilon ? intensity_refl : 0.0f, r.frequency, r.distance_traveled, 0 };
    // writeTextToFile(index_num, r.direction.safeNorm(), "direction_norm2.txt");
    // writeTextToFile(index_num, random_normal.safeNorm(), "random_norm2.txt");
    // writeTextToFile(index_num, refraction_angle, "refraction_angle_snell2.txt");
    // writeTextToFile(index_num, incidence_angle, "reflection_angle_snell2.txt");
    // index_num++;
    return { back_to_transducer_intensity, reflection_ray, refraction_ray };
    
}

units::length::millimeter_t ray_physics::segment_length_in_mm(const btVector3 & v1, const btVector3 & v2)
{
    using namespace std;

    auto x_dist = abs(v1.getX() - v2.getX());
    auto y_dist = abs(v1.getY() - v2.getY());
    auto z_dist = abs(v1.getZ() - v2.getZ());

    //mm
    return units::length::millimeter_t(sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2)));

}
void ray_physics::travel(ray & r, units::length::millimeter_t mm)
{
    r.distance_traveled = r.distance_traveled + mm;
    r.intensity = r.intensity * std::exp(-r.media.attenuation*(mm.to<float>()*0.1f)*r.frequency * 0.1/* dB = 10 log()*/); 
}

bool ray_physics::should_travel(const ray & r)
{
    return r.depth < r.max_depth;
}

float ray_physics::max_ray_length(const ray & r)
{
    // Equation (2), following Beer-Lambert law
    return 10.f /*<- cm to mm*/ * 10.f /* dB = 10*log()*/ *std::log(ray::intensity_epsilon/r.intensity) / (-r.media.attenuation * r.frequency);
}

btVector3 ray_physics::snells_law(const btVector3 & ray_direction, const btVector3 & surface_normal, float incidence_angle, float refraction_angle, float refr_ratio)
{
    // For more details, read https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
    const btVector3 & l = ray_direction;
    const btVector3 & n = surface_normal;
    const float c = incidence_angle;
    const float r = refr_ratio;
    btVector3 temp = ( r * l + (r*c - refraction_angle) * n );

    return temp.safeNormalize();
    
}

float ray_physics::reflection_intensity(const float intensity_in, const float media_1, const float incidence_angle, const float media_2, const float refracted_angle)
{
    const auto && num = media_1 * incidence_angle - media_2 * refracted_angle;
    const auto && denom = media_1 * incidence_angle + media_2 * refracted_angle;

    return intensity_in * pow(num/denom, 2);
}

float ray_physics::reflected_intensity(const float ray_intensity, const float incidence_angle, const material & ray_media, const material & colliding_media)
{
    // Eq. 10 in Burger13
    constexpr auto small_reflections_enhancement_factor = 0.1;
    // TODO: adjust parameter
    constexpr auto custom_reflection_enhancement_factor = 0.2; // we made this up

    const auto specular_factor = std::pow(incidence_angle, colliding_media.specularity);
    const auto impedance_factor = std::pow(( (colliding_media.impedance - ray_media.impedance)
                                            /(colliding_media.impedance + ray_media.impedance)),2);
    const auto intensity = std::pow(ray_intensity, small_reflections_enhancement_factor);

    //std::cout << media_1.impedance << " " << media_2.impedance << std::endl;
    //std::cout << ray_intensity << ", " << specular_factor << " * " << impedance_factor << " * " << intensity << " = " << std::abs(specular_factor * impedance_factor * intensity) << std::endl;
    //return std::pow(std::abs(specular_factor * impedance_factor * ray_intensity), small_reflections_enhancement_factor);
    return std::abs(specular_factor * std::pow(impedance_factor, custom_reflection_enhancement_factor) * intensity);
}
// describes the determination of a random unit vector around v with given polar angle theta.
btVector3 ray_physics::random_unit_vector(btVector3 v, float cos_theta)
{
    bool flag = false;
    float px, py,p;
    do
    {
        //Generating a random point within a circle (uniformly)
        //
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        double a = distribution(generator) * 2 * M_PI;
        double r = 0.5 * sqrt(distribution(generator));
        //in Cartesian coordinates
        px = r * cos(a);
        py = r * sin(a);
        p = px * px + py * py;
    } while (! (p <= 0.25) );
    
    float vx = v.getX();
    float vy = v.getY();
    float vz = v.getZ();
    if ( abs(vx) > abs(vy) )
    {
        vx = vy;
        vy = v.getX();
        flag = true;
    }
    float b = 1 - vx * vx;
    float radicando = 1 - cos_theta * cos_theta;
    radicando = radicando / (p * b);
    float c = sqrt(radicando);
    px = px * c;
    py = py * c;
    float d = cos_theta - vx * px;
    float wx = vx * cos_theta - b * px;
    float wy = vy * d + vz * py;
    float wz = vz * d - vy * py;
    if (flag)
    {
        float aux = wy;
        wy = wx;
        wx = aux;
    }
    return btVector3 (wx,wy,wz);
}

float ray_physics::power_cosine_variate(int v)
{
    //std::default_random_engine generator;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    // std::cout << "aleatorio : " << number << std::endl;
    int indice = v + 1;
    float exponente = (double)1.0 / indice; //TODO: the diffuse or specular should be related to the frequency; cos(theta) should depend on the incidence angle???
    return pow(number, exponente);
}

void ray_physics::writeTextToFile(float index, float text, const char* filename)
    {
    std::ofstream fout(filename, std::ios::app);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }
    if(std::to_string(text) == " ")
    {
        std::cout << "right here" << std::endl;
    }
     fout<<index <<"\t" << std::to_string(text) <<"\t";
    
    fout.close();
    }