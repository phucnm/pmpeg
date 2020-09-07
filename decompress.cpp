//
//  pjpeg_decompress.cpp
//  Toy Image Compressor
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include "input_stream.hpp"
#include "bitmap_image.hpp"
#include "common.hpp"
#include "dct.hpp"
#include "yuv_stream.hpp"
#include "motion_compensation.hpp"
#include <chrono>

const int PLANE_Y = 0;
const int PLANE_CB = 1;
const int PLANE_CR = 2;

//read next 8x8 quantized block
inline void read_block(InputBitStream &input_stream, vector<vector<int>>& block) {
    //Experiment with huffman decoding
    //First read number of symbols
    int n_syms = (int)input_stream.read_byte();
    unordered_map<int, u32> lengths_map;
    for (int i = 0; i < n_syms; i++) {
        int symbol = (int)input_stream.read_u32();
        u32 length = (u32)input_stream.read_byte();
        lengths_map[symbol] = length;
    }

    auto code_table = code_lengths_to_code_table(lengths_map);
    // code: symbol
    unordered_map<pair<u16, u16>, int, hash_pair> lut;
    for (auto c: code_table) {
        lut[c.second] = c.first;
    }
    //Start reading bitstream
//    int count = 64;
    int decoded_size = 0;
    for (unsigned int y = 0; y < 8; y++)
        for (unsigned int x = 0; x < 8; x++) {
            //Read 64 coefs
            int code = 0;
            bool found = false;
            int len = 0;
            while (!found) {
                int read = input_stream.read_bit();
                code = code << 1 | read;
                len++;
                auto pair = make_pair(len, code);
                if (lut.find(pair) != lut.end()) {
                    found = true;
                    block.at(y).at(x) = lut[pair];
                    decoded_size += len;
                }
            }
        }
}

inline void transform_and_fill_zzblock(YUVFrame420& frame,
                                       YUVFrame420& prev_frame,
                                       const vector<vector<int>>& quantize_vector,
                                       const vector<int>& zz_block,
                                       int block_h_idx, int block_w_idx,
                                       int plane_type,
                                       bool is_intra,
                                       Vector2D motion_vector = Vector2D{0, 0}) {
    auto block = create_2d_vector<int>(8, 8);
    zigzag_inflate(zz_block, block);
    auto q_output = create_2d_vector<double>(8, 8);
    dequantize(block, quantize_vector, q_output);
    dct8_2d_inverse_transform(q_output);
    if (plane_type == PLANE_Y) {
        //put idct block to the plane matrix
        for (unsigned int h = 0; h < 8; h++)
            for (unsigned int w = 0; w < 8; w++) {
                u32 x = block_w_idx * 8 + w;
                u32 y = block_h_idx * 8 + h;
                if (x < frame.get_width() && y < frame.get_height()) {
                    unsigned char ref_val = is_intra ? 0 : prev_frame.Y(motion_vector.x, motion_vector.y);
                    frame.Y(x, y) = round_and_clamp_to_char(q_output.at(h).at(w) + ref_val);
                }
            }
    } else if (plane_type == PLANE_CB) {
        //put idct block to the plane matrix
        for (unsigned int h = 0; h < 8; h++)
            for (unsigned int w = 0; w < 8; w++) {
                u32 x = block_w_idx * 8 + w;
                u32 y = block_h_idx * 8 + h;
                if (x < frame.get_width() / 2 && y < frame.get_height() / 2) {
                    unsigned char ref_val = is_intra ? 0 : prev_frame.Cb(motion_vector.x/2, motion_vector.y/2);
                    frame.Cb(x, y) = round_and_clamp_to_char(q_output.at(h).at(w) + ref_val);
                }
            }
    } else if (plane_type == PLANE_CR) {
        //put idct block to the plane matrix
        for (unsigned int h = 0; h < 8; h++)
            for (unsigned int w = 0; w < 8; w++) {
                u32 x = block_w_idx * 8 + w;
                u32 y = block_h_idx * 8 + h;
                if (x < frame.get_width() / 2 && y < frame.get_height() / 2) {
                    unsigned char ref_val = is_intra ? 0 : prev_frame.Cr(motion_vector.x/2, motion_vector.y/2);
                    frame.Cr(x, y) = round_and_clamp_to_char(q_output.at(h).at(w) + ref_val);
                }
            }
    }
}

inline void read_plane2(YUVFrame420 &frame,
                        YUVFrame420& prev_frame,
                        InputBitStream& input_stream,
                        u32 width, u32 height,
                        const vector<vector<int>>& quantize_vector,
                        int plane_type,
                        bool is_intra,
                        const vector<vector<Vector2D>> motion_vects = vector<vector<Vector2D>>(),
                        int plane_scale = 1) {
    int n_syms = (int)input_stream.read_byte();
    unordered_map<RunSize, u32> lengths_map;
    for (unsigned int i = 0; i < n_syms; i++) {
        int run = input_stream.read_byte();
        int size = input_stream.read_byte();
        auto symbol = RunSize{run, size};
        u32 length = (u32)input_stream.read_byte();
        lengths_map[symbol] = length;
    }

    auto code_table = code_lengths_to_code_table(lengths_map);
    // code: symbol
    unordered_map<pair<u16, u16>, RunSize, hash_pair> lut;
    for (auto c: code_table) {
        lut[c.second] = c.first;
    }

    for (unsigned int bh = 0; bh < height / 8; bh++)
        for (unsigned int bw = 0; bw < width / 8; bw++) {
            Vector2D motion_vector = Vector2D{0, 0};

            //What's this macro block index?
            if (!is_intra) {
                int mb_x = bw * plane_scale / 2;
                int mb_y = bh * plane_scale / 2;
                if (mb_y < motion_vects.size() && mb_x < motion_vects[0].size()) {
                    Vector2D mv = motion_vects.at(mb_y).at(mb_x);
                    if (mv.skip) {
                        //skip this block
                        continue;
                    } else {
                        motion_vector = mv;
                    }
                } else {
                    continue;
                }
            }
            //after zigzag ordering
            vector<int> zigzag_block(64, 0);
            int num_decoded = 0;
            int count = 0;
            bool eob = false;
            while (!eob) {
                //Read block
                int code = 0;
                bool found = false;
                int len = 0;
                while (!found) {
                    int read = input_stream.read_bit();
                    code = code << 1 | read;
                    len++;
                    auto pair = make_pair(len, code);
                    if (lut.find(pair) != lut.end()) {
                        found = true;
                        RunSize rs = lut[pair];
                        num_decoded += rs.run + 1;
                        if ((rs.run == 0 && rs.size == 0)) {
                            eob = true;
                            transform_and_fill_zzblock(frame, prev_frame, quantize_vector, zigzag_block, bh, bw, plane_type, is_intra, motion_vector);
                            break;
                        }
                        //read the diff value
                        int val = input_stream.read_bits(rs.size);
                        int low_positive_bound = 1<<(rs.size - 1);
                        if (val < low_positive_bound) {
                            int min_val_in_bits_length = -((1<<rs.size) - 1);
                            val = min_val_in_bits_length + val;
                        }

                        int run = rs.run;
                        while (--run >= 0) {
                            zigzag_block[count++] = 0;
                        }
                        zigzag_block[count++] = val;
                        if (num_decoded == 64) {
                            eob = true;
                            transform_and_fill_zzblock(frame, prev_frame, quantize_vector, zigzag_block, bh, bw, plane_type, is_intra, motion_vector);
                            break;
                        }
                    }
                }
            }

        }
}

inline void decompress(YUVFrame420& frame,
                       YUVFrame420& prev_frame,
                       InputBitStream& input_stream,
                       u32 width, u32 height,
                       const vector<vector<int>>& quantize_vector,
                       bool is_intra,
                       const vector<vector<Vector2D>> motion_vects = vector<vector<Vector2D>>()) {
    int padded_height = next_mul(height, 8);
    int padded_width = next_mul(width, 8);

    int padded_c_h = next_mul((height + 1) / 2, 8);
    int padded_c_w = next_mul((width + 1) / 2, 8);

//    auto Y = create_2d_vector<int>(padded_height,padded_width);
//    auto Y_plane = create_2d_vector<unsigned char>(padded_height, padded_width);
//    auto Cb_plane = create_2d_vector<unsigned char>(padded_c_h, padded_c_w);
//    auto Cr_plane = create_2d_vector<unsigned char>(padded_c_h, padded_c_w);
    read_plane2(frame, prev_frame, input_stream, padded_width, padded_height, quantize_vector, PLANE_Y, is_intra, motion_vects, 1);
    read_plane2(frame, prev_frame, input_stream, padded_c_w, padded_c_h, quantize_vector, PLANE_CB, is_intra, motion_vects, 2);
    read_plane2(frame, prev_frame, input_stream, padded_c_w, padded_c_h, quantize_vector, PLANE_CR, is_intra, motion_vects, 2);
    input_stream.flush_to_byte();
}

inline void decompress_P_frame(YUVFrame420& frame, YUVFrame420& prev_frame, InputBitStream& input_stream, u32 width, u32 height, const vector<vector<int>>& quantize_vector) {
    auto mo_vects = create_2d_vector<Vector2D>(height / 16, width / 16);
    for (unsigned int y = 0; y < height/16; y++) {
        for (unsigned int x = 0; x < width / 16; x++) {
            int has_motion = input_stream.read_bit();
            if (has_motion) {
                u16 mo_x = input_stream.read_u16();
                u16 mo_y = input_stream.read_u16();
                mo_vects.at(y).at(x) = Vector2D{mo_x, mo_y};
            } else {
                mo_vects.at(y).at(x) = Vector2D{0, 0, .skip = true};
            }
        }
    }
    decompress(frame, prev_frame, input_stream, width, height, quantize_vector, false, mo_vects);
}

void do_decompress() {
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
}
