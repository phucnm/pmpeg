/* uvid_decompress.cpp
   CSC 485B/578B/SENG 480B - Data Compression - Summer 2020

   Starter code for Assignment 5
   
   This placeholder code reads the (basically uncompressed) data produced by
   the uvid_compress starter code and outputs it in the uncompressed 
   YCbCr (YUV) format used for the sample video input files. To play the 
   the decompressed data stream directly, you can pipe the output of this
   program to the ffplay program, with a command like 

     ffplay -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 - 2>/dev/null
   (where the resolution is explicitly given as an argument to ffplay).

   B. Bird - 07/15/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include <tuple>
#include "input_stream.hpp"
#include "yuv_stream.hpp"
#include "decompress.cpp"

const int GOP_COUNT = 3;

int main(int argc, char** argv){

    //Note: This program must not take any command line arguments. (Anything
    //      it needs to know about the data must be encoded into the bitstream)
    
    InputBitStream input_stream {std::cin};

    u32 height {input_stream.read_u32()};
    u32 width {input_stream.read_u32()};

    YUVStreamWriter writer {std::cout, width, height};

    //Read quality code
    vector<vector<int>> quantize_vector = const_medium_quantize_vector;
    int quality_code = input_stream.read_bits(2);
    float quality_factor = 1;
    if (quality_code == LOW_QUALITY_CODE) {
        quality_factor = 2;
    } else if (quality_code == HIGH_QUALITY_CODE) {
        quality_factor = 0.5;
    }
    for (int x = 0; x < 8; x++)
        for (int y = 0; y < 8; y++)
            quantize_vector.at(x).at(y) *= quality_factor;

    int count = 0;
    YUVFrame420 f{0, 0};
    YUVFrame420& last_decompressed_frame = f;

    while (input_stream.read_byte()){
        YUVFrame420& frame = writer.frame();
        if (count % GOP_COUNT == 0) {
            decompress(frame, last_decompressed_frame, input_stream, width, height, quantize_vector, true);
        } else {
            frame.copy_from_frame(last_decompressed_frame);
            decompress_P_frame(frame, last_decompressed_frame, input_stream, width, height, quantize_vector);
        }

        count++;
        writer.write_frame();
        last_decompressed_frame = frame;
    }
    input_stream.flush_to_byte();

    return 0;
}
