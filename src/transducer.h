#ifndef TRANSDUCER_H
#define TRANSDUCER_H

#include <units/units.h>
#include <LinearMath/btVector3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

//#define M_PI 3.14159

/**
 * Objective: define a class transducer, which could generate the position and direction of each ray (sound wave)
 * 
 * Note: pay more attention to how to assign the angles and position of transducer in 3D model scene.
 * 
 * 
 */


template<size_t transducer_elements, size_t sample_count>
// template<size_t transducer_elements>
class transducer
{
public:
    struct transducer_element
    {
        btVector3 position; // unit: mm. 
        btVector3 direction; // unit
    };

    transducer(const float frequency/*MHz*/, const units::length::centimeter_t radius, units::length::millimeter_t transducer_element_separation,
               const btVector3 & position, const std::array<units::angle::degree_t, 3> & angles) :
        frequency(frequency),
        radius(radius),
        position(position), //unit: mm. position of the transducer with respect to world coordinate system (3D scene)
        angles(angles) // unit: degree. angles of the transducer with respect ot world coordinate system (3D scene), rotation order Z-> X -> Y
    {
        using namespace units::angle;
        using namespace units::literals;

        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<double> distribution(0.0, 0.024);

        assert(transducer_element_separation * transducer_elements < M_PI * radius);

        radian_t x_angle { angles[0] };
        radian_t y_angle { angles[1] };
        radian_t z_angle { angles[2] };

        auto amp = transducer_element_separation / radius;
        const radian_t amplitude { amp.to<float>() }; // angle covered by a single TE
        const radian_t angle_center_of_element { amplitude / 2.0f };

        radian_t angle = -(amplitude * transducer_elements / 2) + angle_center_of_element;

        for (size_t t = 0; t < transducer_elements; t++)
        {
            auto & samples_elements = elements[t];

            for(size_t i = 0; i < sample_count; i++)
            {
                double delta_x = distribution(generator);
                // double delta_x = 0;
                double delta_z = distribution(generator);
                // double delta_z = 0;
                samples_elements[i] = transducer_element
                {   
                    position + btVector3 (delta_x, 10*radius.to<float>(), delta_z).rotate(btVector3(0,0,1), -angle.to<float>()).rotate(btVector3(0,0,1), z_angle.to<float>())
                                                                                                                                        .rotate(btVector3(1,0,0), x_angle.to<float>())
                                                                                                                                        .rotate(btVector3(0,1,0), y_angle.to<float>()), // position
                    btVector3 ( std::sin(angle.to<float>()), std::cos(angle.to<float>()), 0 ).rotate(btVector3(0,0,1), z_angle.to<float>())
                                                                                         .rotate(btVector3(1,0,0), x_angle.to<float>())
                                                                                         .rotate(btVector3(0,1,0), y_angle.to<float>())  // direction

                    // position + 10*radius.to<float>() * btVector3 ( std::sin(angle.to<float>()), std::cos(angle.to<float>()), 0 ).rotate(btVector3(0,0,1), z_angle.to<float>())
                    //                                                                                                      .rotate(btVector3(1,0,0), x_angle.to<float>())
                    //                                                                                                      .rotate(btVector3(0,1,0), y_angle.to<float>()), // position
                    // btVector3 ( std::sin(angle.to<float>()), std::cos(angle.to<float>()), 0 ).rotate(btVector3(0,0,1), z_angle.to<float>())
                    //                                                                      .rotate(btVector3(1,0,0), x_angle.to<float>())
                    //                                                                      .rotate(btVector3(0,1,0), y_angle.to<float>())  // direction
                };
     

            }
            angle = angle + amplitude;
            // std::cout << angle << std::endl;
            
        }

    }

    transducer_element element(size_t index_ray, size_t index_sample) const
    {   
        auto samples_element = elements[index_ray];
        return samples_element[index_sample];
        // return elements[index_ray][index_sample];
        // return elements.at(i);
    }

    void print(bool direction) const
    {
        auto print_vec = [](const auto & v)
        {
            std::cout << v.x() << "," << v.z() << std::endl;
        };

        for (auto & element : elements)
        {
            print_vec(direction? element.direction : element.position);
        }
    }

    const float frequency;

    const btVector3 position, direction;
    const std::array<units::angle::degree_t, 3> angles;

private:
    const units::length::millimeter_t radius;

    // std::array<transducer_element, transducer_elements> elements;
    std::array<std::array<transducer_element, sample_count>, transducer_elements> elements;

};

#endif // TRANSDUCER_H
