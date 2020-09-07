/* uvid_compress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5

   Reads video data from stdin in uncompresed YCbCr (YUV) format 
   (With 4:2:0 subsampling). To produce this format from 
   arbitrary video data in a popular format, use the ffmpeg
   tool and a command like 

     ffmpeg -i videofile.mp4 -f rawvideo -pixel_format yuv420p - 2>/dev/null | ./this_program <width> height>

   Note that since the width/height of each frame is not encoded into the raw
   video stream, those values must be provided to the program as arguments.

   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "output_stream.hpp"
#include "yuv_stream.hpp"
#include "compress.cpp"
#include "dct.hpp"
#include <chrono>
#include <random>

//static std::default_random_engine randGen((std::random_device())());



int main(int argc, char** argv){

    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <width> <height> <low/medium/high>" << std::endl;
        return 1;
    }
    u32 width = std::stoi(argv[1]);
    u32 height = std::stoi(argv[2]);
    std::string quality{argv[3]};

    YUVStreamReader reader {std::cin, width, height};
    OutputBitStream output_stream {std::cout};

    output_stream.push_u32(height);
    output_stream.push_u32(width);

    //Prepare quantization vector
    // For low just /2, high * 2 otherwise leave it as is.
    vector<vector<int>> quantize_vector = const_medium_quantize_vector;
    int quality_code = MEDIUM_QUALITY_CODE;
    float quality_factor = 1;
    if (quality.compare("low") == 0) {
        quality_factor = 2;
        quality_code = LOW_QUALITY_CODE;
    } else if (quality.compare("high") == 0) {
        quality_factor = 0.5;
        quality_code = HIGH_QUALITY_CODE;
    }
    for (unsigned int x = 0; x < 8; x++)
        for (unsigned int y = 0; y < 8; y++)
            quantize_vector.at(x).at(y) *= quality_factor;
    //Output the quality into bitstream
    output_stream.push_bits(quality_code, 2);
    YUVFrame420 f{0, 0};
    YUVFrame420& last_decompressed_frame = f;

    int count = 0;
    while (reader.read_next_frame()){
        output_stream.push_byte(1); //Use a one byte flag to indicate whether there is a frame here
        YUVFrame420& frame = reader.frame();
        if (count % GOP_COUNT == 0) {
            compress_frame(frame, last_decompressed_frame, width, height, output_stream, quantize_vector, true);
        } else {
            compress_P_frames(frame, last_decompressed_frame, width, height, output_stream, quantize_vector);
        }
        ++count;
        last_decompressed_frame = frame;
    }

    output_stream.push_byte(0); //Flag to indicate end of data
    output_stream.flush_to_byte();

    return 0;
}
