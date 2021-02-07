#ifndef RFIMAGE_H
#define RFIMAGE_H

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <array>
#include <units/units.h>
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex> 

#include "psf.h"
using namespace std;

//macros for real and imaginary parts
#define REAL 0
#define IMAG 1

/**
 * Radio-frequency image.
 *
 * Stores the resulting echoes of the ultrasound wave interacting with tissues.
 * Each column stands for information from a single transducer element.
 * The number of rows is calculated automatically according to maximum travel
 * time of the ultrasound pulse, psf's axial resolution, and average speed of sound.
 */
// columns: the number of transducer elements
// max_rows: 


template <unsigned int columns, unsigned int max_travel_time /*μs*/, unsigned int axial_resolution /*μm*/, unsigned int speed_of_sound = 1500 /*μm/μs*/>
class rf_image
{
public:
    rf_image(units::length::millimeter_t radius, units::angle::radian_t angle) :
        intensities(max_rows, columns, CV_32FC1),
        conv_axial_buffer(max_rows, columns, CV_32FC1),
        scan_converted(700,900, CV_32FC1)
    {
        std::cout << "rf_image: " << max_rows << ", " << columns << std::endl;
        create_mapping(radius, angle, columns, max_rows);
    }
    
    void add_echo(const unsigned int column, const float echo, const units::time::microsecond_t micros_from_source) // micros_from_source should be the time from the starting of US image to the point
    {
        const units::dimensionless::dimensionless_t row = micros_from_source / (axial_resolution_ / speed_of_sound_);

        if (row < max_rows)
        {
            float echo_previous = intensities.at<float>(row, column);
            if (echo_previous < echo) // chose the max one as the pixel value
            {
                intensities.at<float>(row, column) = echo; 
            } 
            // intensities.at<float>(row, column) += echo;
        }
    }

    // get the delta time that represents a pixel (row resolution in time)
    constexpr units::time::microsecond_t get_dt() const
    {
        return axial_resolution / speed_of_sound_;
    }

    constexpr units::time::microsecond_t micros_traveled(units::length::micrometer_t microm_from_source) const
    {
        return microm_from_source / speed_of_sound_;
    }

    // Transforms the rf image by doing a fast approximation of the envelope function
    void envelope()
    {
        // Travel through each column looking for concave peaks.
        // Then, recalculate the points between the peaks as the linear interpolation of the absolute values of those peaks.
        // This should work as a fast approximation of the hilbert transform over the rf signal.

        // for (size_t column = 0; column < columns; column++)
        // {
        //     bool ascending = intensities.at<float>(0, column) < intensities.at<float>(1, column);
        //     size_t last_peak_pos = 0;
        //     float last_peak = intensities.at<float>(last_peak_pos, column);
        //     for (size_t i = 1; i < max_rows-1; i++)
        //     {
        //         if (intensities.at<float>(i, column) < intensities.at<float>(i+1, column))
        //         {
        //             ascending = true;
        //         }
        //         else if (ascending)
        //         // if it was ascending and now descended, we found a concave point at i
        //         {
        //             ascending = false;
        //             const float new_peak = std::abs(intensities.at<float>(i, column));

        //             // lerp last_peak -> new_peak over last_peak_pos -> i (new_peak_pos)
        //             for (size_t j = last_peak_pos; j < i; j++)
        //             {
        //                 const float alpha = (static_cast<float>(j) - static_cast<float>(last_peak_pos)) /
        //                                     (static_cast<float>(i) - static_cast<float>(last_peak_pos));

        //                 intensities.at<float>(j, column) = last_peak * (1-alpha) + new_peak * alpha;
        //             }

        //             last_peak_pos = i;
        //             last_peak = new_peak;
        //         }
        //     }
        // }
        int N = intensities.rows;
        fftw_complex y[N];
        
        for (int i = 0; i < intensities.cols; i++)
        {
            cv::Mat in;
            intensities.col(i).copyTo(in);
            cv::Scalar meanValue = cv::mean(in);
            in = in - meanValue;
            hilbert(in, y);
            for (int j = 0; j < N; j++)
            {
                complex<double> mycomplex(y[j][REAL], y[j][IMAG]);
                intensities.at<float>(j,i) = meanValue.val[0] + abs(mycomplex);
            }
        }
        
    }

    template <typename psf_>
    void convolve(const psf_ & p)
    {
        // Convolve using only axial kernel and store in intermediate buffer
        int half_axial = p.get_axial_size()/2;
        for (int col = 0; col < intensities.cols; col++) //each column is a different TE
        {
            for (int row = half_axial; row < intensities.rows - half_axial; row++) // each row holds information from along a ray
            {
                float convolution = 0;
                for (int kernel_i = 0; kernel_i < p.get_axial_size(); kernel_i++)
                {
                    convolution += intensities.at<float>(row - half_axial+ kernel_i, col) * p.axial_kernel[kernel_i];
                }
                conv_axial_buffer.at<float>(row,col) = convolution;
            }
        }

        // Convolve intermediate buffer using lateral kernel
        int half_lateral =  p.get_lateral_size() / 2;
        for (int row = 0; row < conv_axial_buffer.rows; row++) // each row holds information from along a ray
        {
            for (int col = half_lateral; col < conv_axial_buffer.cols - half_lateral; col++) //each column is a different TE
            {
                float convolution = 0;
                for (int kernel_i = 0; kernel_i < p.get_lateral_size(); kernel_i++)
                {
                    convolution += conv_axial_buffer.at<float>(row, col - half_lateral + kernel_i) * p.lateral_kernel[kernel_i];
                }
                // scale image [-6 6] map to [0, 0.8]

                if (convolution < -0.06)
                {
                    intensities.at<float>(row,col) = 0;
                }
                else
                {
                    intensities.at<float>(row, col) = (convolution + 0.06)/1.06;
                }        
            }
        }

    }

    void postprocess()
    {
        

        // cv::imwrite("prelog_rf.png", intensities);
        writeMatToFile(intensities, "prelog_rf.txt");
        cv::Scalar temp = cv::mean(intensities);
        double meanValue = temp.val[0];
        for (size_t i = 0; i < max_rows * columns; i++)
        {
            
            intensities.at<float>(i) = 20*std::log10(intensities.at<float>(i)+0.1 / (meanValue + 0.1));
            // intensities.at<float>(i) = std::log10(intensities.at<float>(i)+1)/std::log10(max+1);

        }
        double min, max;
        cv::minMaxLoc(intensities, &min, &max);
        for (size_t i = 0; i < max_rows * columns; i++)
        {   
            float temp = intensities.at<float>(i);
            if (temp <= 0)
            {
                intensities.at<float>(i) = temp*(-0.5)/min + 0.5; 
            }
            else
            {
                intensities.at<float>(i) = temp * 0.5/max + 0.5;
            }
        }
        double contrast = 0.4; // map [0,1] to [0,0.5] to [0,1]
         for (size_t i = 0; i < max_rows * columns; i++)
        {   
            float temp = intensities.at<float>(i);
            if (temp <= contrast)
            {
                intensities.at<float>(i) = (temp - 0.5)/0.5 + 1; 
            }
            else
            {
                intensities.at<float>(i) = 1;
            }
        }
        
        writeMatToFile(intensities, "postlog_rf.txt");
        //cv::imwrite("postlog_rf.png", intensities);
        // apply scan conversion using preprocessed mapping
        constexpr float invalid_color = 0.0f;
        cv::remap(intensities, scan_converted, map_y, map_x, CV_INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(invalid_color)); 
        writeMatToFile(scan_converted, "Simulated_US.txt");
    }

    void save(const std::string & filename) const
    {
        cv::imwrite(filename, intensities);
    }

    void show() const
    {
        //cv::namedWindow("Display window", cv::WINDOW_AUTOSIZE );
        //cv::imshow("Display window", image );

        cv::namedWindow("Scan Converted", cv::WINDOW_AUTOSIZE );
        cv::imshow("Scan Converted", scan_converted );

        cv::waitKey(16);
    }

    void clear()
    {
        intensities.setTo(-2.0f); // Tricky point
    }

    void print(size_t column) const
    {
        for (size_t i = 0; i < max_rows; i++)
        {
            std::cout << intensities.at<float>(i, column) << ", ";
        }
        std::cout << std::endl;
    }


    void writeMatToFile(cv::Mat& m, const char* filename)
    {
    std::ofstream fout(filename);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }

    for(int i=0; i<m.rows; i++)
    {
        for(int j=0; j<m.cols; j++)
        {
            fout<<m.at<float>(i,j)<<"\t";
        }
        fout<<endl;
    }

    fout.close();
    }

    void hilbert(cv::Mat in, fftw_complex* out)
    {
        int N = intensities.rows;
        // copy the data to the complex array
        for (int i = 0; i < N; ++i) {
            out[i][REAL] = in.at<float>(i,0);
            out[i][IMAG] = 0;
        }
        // creat a DFT plan and execute it
        fftw_plan plan = fftw_plan_dft_1d(N, out, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        // destroy a plan to prevent memory leak
        fftw_destroy_plan(plan);
        int hN = N >> 1; // half of the length (N/2)
        int numRem = hN; // the number of remaining elements
        // multiply the appropriate value by 2
        //(those should multiplied by 1 are left intact because they wouldn't change)
        for (int i = 1; i < hN; ++i) {
            out[i][REAL] *= 2;
            out[i][IMAG] *= 2;
        }
        // if the length is even, the number of the remaining elements decrease by 1
        if (N % 2 == 0)
            numRem--;
        else if (N > 1) {
            out[hN][REAL] *= 2;
            out[hN][IMAG] *= 2;
        }
        // set the remaining value to 0
        // (multiplying by 0 gives 0, so we don't care about the multiplicands)
        memset(&out[hN + 1][REAL], 0, numRem * sizeof(fftw_complex));
        // creat a IDFT plan and execute it
        plan = fftw_plan_dft_1d(N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        // do some cleaning
        fftw_destroy_plan(plan);
        fftw_cleanup();
        // scale the IDFT output
        for (int i = 0; i < N; ++i) {
            out[i][REAL] /= N;
            out[i][IMAG] /= N;
        }

    }   


private:
    // Create constexpr variables to use type operations later
    static constexpr units::velocity::meters_per_second_t speed_of_sound_ = units::velocity::meters_per_second_t(speed_of_sound);
    static constexpr units::length::micrometer_t axial_resolution_ = units::length::micrometer_t(axial_resolution);

    static constexpr unsigned int max_rows = (speed_of_sound * max_travel_time) / axial_resolution;

    // fills map_x and map_y
    void create_mapping(units::length::millimeter_t radius, units::angle::radian_t total_angle, unsigned int rf_width, unsigned int rf_height)
    {
        // ratio to convert from mm to px
        float ratio = (max_travel_time * speed_of_sound * 0.001f + radius.to<float>() - radius.to<float>() * std::cos(total_angle.to<float>()/2.0)) / scan_converted.rows;

        // distance to transducer center going from the edge of the scan_converted image
        units::length::millimeter_t shift_y = radius * std::cos(total_angle.to<float>() / 2.0f);

        // horizontal center of the image
        float half_width = (float)scan_converted.cols / 2.0f;

        map_x.create(scan_converted.size(), CV_32FC1);
        map_y.create(scan_converted.size(), CV_32FC1);

        for( int j = 0; j < scan_converted.cols; j++ )
        {
            for( int i = 0; i < scan_converted.rows; i++ )
            {
                float fi = static_cast<float>(i)+shift_y.to<float>()/ratio;
                float fj = static_cast<float>(j)-half_width;

                // distance from the point i,j to the center of the transducer, located in -shift_y,half_width
                float r = std::sqrt(std::pow(fi,2.0f) + std::pow(fj,2.0f));

                // angle of the vector (i+shift_y, j-half_width) against the half center line
                units::angle::radian_t angle = units::angle::radian_t(std::atan2(fj, fi));

                // Invalid values are OK here. cv::remap checks for them and assigns them a special color.
                // Note: in scan converted image, given that the pixel size is "ratio".
                map_x.at<float>(i,j) = (r*ratio-radius.to<float>())/(max_travel_time*speed_of_sound*0.001f) * (float)rf_height;
                map_y.at<float>(i,j) = ((angle - (-total_angle/2)) / (total_angle)) * (float)rf_width;
            }
        }
    }

    cv::Mat intensities;
    cv::Mat scan_converted;// note: scan_converted is the image we usually see from US machine
    cv::Mat conv_axial_buffer;

    // Mapping used for scan conversion
    cv::Mat map_x, map_y;
};

// http://stackoverflow.com/a/22414046
template <unsigned int columns, unsigned int max_travel_time, unsigned int axial_resolution, unsigned int speed_of_sound>
constexpr units::velocity::meters_per_second_t rf_image<columns, max_travel_time, axial_resolution, speed_of_sound>::speed_of_sound_;

template <unsigned int columns, unsigned int max_travel_time, unsigned int axial_resolution, unsigned int speed_of_sound>
constexpr units::length::micrometer_t rf_image<columns, max_travel_time, axial_resolution, speed_of_sound>::axial_resolution_;

#endif // RFIMAGE_H
