#include "transducer.h"
#include "scene.h"
#include "volume.h"
#include "psf.h"
#include "rfimage.h"


#include <cmath>
#include <iostream>
#include <fstream>
#include <units/units.h>
#include <nlohmann/json.hpp>

using namespace units::literals;
using namespace units::velocity;
using namespace units::length;
using namespace units::time;
using namespace units::angle;

constexpr meters_per_second_t speed_of_sound = 1540_m / 1_s; // [μm/μs], [m/s]
constexpr float transducer_frequency = 4.5f; // [Mhz]
// TODO: which way to get axial resolution is more persuasive?
constexpr millimeter_t axial_resolution = millimeter_t(1.45f / transducer_frequency); // [mm], the division can be deduced from Burger13 
constexpr size_t transducer_elements = 1024;
constexpr size_t samples_count = 20; 
constexpr radian_t transducer_amplitude = 70_deg;
constexpr centimeter_t transducer_radius = 3_cm;
constexpr centimeter_t ultrasound_depth = 15_cm; // [15cm -> μm]

constexpr microsecond_t max_travel_time = microsecond_t(ultrasound_depth / speed_of_sound); // [μs]

constexpr unsigned int resolution_axial = 145;//145; // [μm], from Burger13
using psf_ = psf<11, 11, 11, 100>;
using volume_ = volume<256, 100>;
using rf_image_ = rf_image<transducer_elements, max_travel_time.to<unsigned int>(), static_cast<unsigned int>(axial_resolution.to<float>()*1000.0f/*mm->μm*/)>;
// using transducer_ = transducer<transducer_elements>;
using transducer_ = transducer<transducer_elements, samples_count>;

std::random_device rd;
std::mt19937 generator(rd());
std::normal_distribution<double> distr(0.0001, 0.01);


int main(int argc, char** argv)

{
    if (argc != 2)
    {
        std::cout << "Incorrect argument list." << std::endl;
        return 0;
    }

    // Step 1: Read 3D mesh model using json
    nlohmann::json json;
    {
        std::ifstream infile { argv[1] };
       // std::ifstream infile { filename };
        json << infile;
    } 

    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << std::endl;
    }

    // Step 2: Initialization (including transducer, scatter distribution, point spread function and radio frequency image)
    // Step 2.1: transducer
    millimeter_t transducer_element_separation = transducer_amplitude.to<float>() * transducer_radius / transducer_elements;
    const auto & t_pos = json.at("transducerPosition"); //unit: mm. position of the transducer with respect to world coordinate system (3D scene)
    const auto & t_dir = json.at("transducerAngles"); // unit: degree. angles of the transducer with respect ot world coordinate system (3D scene), rotation order Z-> X -> Y
    std::array<units::angle::degree_t, 3> transducer_angles = {degree_t((float)t_dir[0]), degree_t((float)t_dir[1]), degree_t((float)t_dir[2])};
    // std::cout << "transducer postition: " << t_pos[0]<<","<< t_pos[1] << "," << t_pos[2]<< std::endl;
    // std::cout << "transducer angles: " << t_dir[0]<<","<< t_dir[1] << "," << t_dir[2]<< std::endl;
    transducer_ transducer(transducer_frequency, transducer_radius, transducer_element_separation,
                           btVector3(t_pos[0], t_pos[1], t_pos[2]), transducer_angles);
    
    // Step 2.2 scatter distribution
    static const volume_ texture_volume;

    // Step 2.3 point spread function
    const float var_x = 0.1; //0.145*0.145
    const float var_y = 0.2; // 0.3*0.3
    const float var_z = 0.2;
    const psf_ psf { transducer_frequency, var_x, var_y, var_z };

    // Step 2.4 radio frequency image
    rf_image_ rf_image { transducer_radius, transducer_amplitude };

    // std::cout << max_travel_time << std::endl;

    try
    {
        scene scene { json, transducer };

        scene.step(1000.0f);
        // while(true)
        //Set transducer movement 
        for(int indexMov = 0; indexMov < 20; indexMov++)
        // int indexMov = 0;
        {
            auto temp_pos = (float)t_pos[2];
            // auto temp_pos = (float)t_pos[2]   - indexMov;
            std::array<units::angle::degree_t, 3> movedAngle = {degree_t((float)t_dir[0]+0), degree_t((float)t_dir[1] + 0), degree_t((float)t_dir[2]+0)};
            transducer_ transducer_temp(transducer_frequency, transducer_radius, transducer_element_separation,
                           btVector3(t_pos[0], t_pos[1], temp_pos), movedAngle);
            rf_image.clear(); // Default set to (-2)

            // auto rays = scene.cast_rays<transducer_elements>();
            auto rays = scene.cast_rays<samples_count, transducer_elements>(transducer_temp);
            // auto rays = scene.cast_rays<samples_count, transducer_elements>(transducer);
            for (unsigned int ray_i = 0; ray_i < rays.size(); ray_i++)
            {
                const auto & ray = rays[ray_i];
                for (unsigned int sample_i = 0; sample_i < samples_count; sample_i++)
                {
                    const auto & sample = ray[sample_i];
                    int num_segment = sample.size();
                    // for (auto & segment : ray)
                    for(auto & segment: sample)
                    {
                        const auto starting_micros = rf_image.micros_traveled(segment.distance_traveled /*mm -> μm*/);
                        const auto distance = scene.distance(segment.from, segment.to); // [mm]
                        const auto steps = distance / axial_resolution;
                        const auto delta_step = axial_resolution.to<float>() * segment.direction;
                        const auto time_step = rf_image.micros_traveled(axial_resolution); // [μs]

                        auto point = segment.from;
                        auto time_elapsed = starting_micros;
                        auto intensity = segment.initial_intensity;

                        for (unsigned int step = 0; step < steps && time_elapsed < max_travel_time; step++)
                        {
                            float scattering = texture_volume.get_scattering(segment.media.mu1, segment.media.mu0, segment.media.sigma, point.x(), point.y(), point.z());
                            // float scattering = 0;
                            // float scatter = intensity * scattering + distr(generator);
                            float scatter = intensity * scattering;
                            // float scatter = scattering;
                            // scatter = 0;


                            rf_image.add_echo(ray_i, scatter, time_elapsed);
                        
                            // Step forward through the segment, decreasing intensity using Beer-Lambert's law
                            point += delta_step;
                            time_elapsed = time_elapsed + time_step;

                            constexpr auto k = 0.1f;
                            intensity *= std::exp(-segment.attenuation * axial_resolution.to<float>()*0.1f * transducer_frequency * k);
                        }

                        // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.
                        // rf_image.add_echo(ray_i, segment.reflected_intensity, starting_micros + time_step * (steps-1));
                        // rf_image.add_echo(ray_i, (segment.reflected_intensity*(-10000)), starting_micros + time_step * (steps-1));
                    }
                    
                }    

            }

            rf_image.convolve(psf);

            for (unsigned int ray_i = 0; ray_i < rays.size(); ray_i++)
            {
                const auto & ray = rays[ray_i];
                for (unsigned int sample_i = 0; sample_i < samples_count; sample_i++)
                {
                    const auto & sample = ray[sample_i];
                    int num_segment = sample.size();
                    // for (auto & segment : ray)
                    for(auto & segment: sample)
                    {
                        const auto starting_micros = rf_image.micros_traveled(segment.distance_traveled /*mm -> μm*/);
                        const auto distance = scene.distance(segment.from, segment.to); // [mm]
                        const auto steps = distance / axial_resolution;
                        const auto time_step = rf_image.micros_traveled(axial_resolution); // [μs]

                        // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.
                        // rf_image.add_echo(ray_i, segment.reflected_intensity, starting_micros + time_step * (steps-1));
                        rf_image.add_echo(ray_i, (segment.reflected_intensity * (1000)), starting_micros + time_step * (steps-1));
                    }
                    
                }    

            }

            rf_image.envelope();
            rf_image.postprocess();
            // rf_image.show();
            
        }
    }
    catch (const std::exception & ex)
    {
        std::cout << "The program found an error and will terminate.\n"
                  << "Reason:\n"
                  << ex.what() << std::endl;
    }

}
