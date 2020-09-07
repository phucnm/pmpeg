//
//  pjpeg_compress.cpp
//
//  Created by Phuc Nguyen on 07/21/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#ifndef COMPRESS_HPP
#define COMPRESS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdint>
#include <algorithm>
#include "output_stream.hpp"
#include "common.hpp"
#include "dct.hpp"
#include <unordered_map>
#include <queue>
#include "yuv_stream.hpp"
#include "motion_compensation.hpp"

/*
 Huffman Tree Node
 */
class Node {
public:
    Node(u32 symbol, u32 value) {
        this->symbol = symbol;
        this->value = value;
    }

    u32 symbol;
    u32 value;
    Node* left;
    Node* right;
};

/*
 Huffman Tree builder.
 Provides a canonical codes builder and package-merge builder as static methods.
 */

class Compare
{
public:
     bool operator ()(const Node* lhs, const Node* rhs) const
       {
           return lhs->value > rhs->value;
       }
};

class HuffmanCodingBuilder {
public:
    //This outputs a map with format [symbol: bit length]
    template<typename T>
    static unordered_map<T, u32> package_merge(const unordered_map<T, int>& freqs, int max_bits) {
        vector<pair<vector<T>, u32>> original_freqs;
        for (auto f: freqs) {
            original_freqs.push_back(make_pair(vector<T>{f.first}, f.second));
        }
        if (original_freqs.empty()) {
            return unordered_map<T, u32>{};
        } else if (original_freqs.size() == 1) {
            auto single_sym = original_freqs[0];
            T sym = single_sym.first[0];
            return unordered_map<T, u32>({ {sym, 1} });
        }
//        vector<pair<u32, u32>> res;
        vector<pair<vector<T>, u32>> pkgs(original_freqs);
        sort(pkgs.begin(), pkgs.end(), [](const pair<vector<T>, u32>& a, const pair<vector<T>, u32>& b) -> bool {
            return a.second < b.second;
        });
        u32 num_pick_items = (u32)(2 * original_freqs.size() - 2);
        while (--max_bits) {
            if (pkgs.size() % 2 == 1) {
                pkgs.pop_back();
            }
            vector<pair<vector<T>, u32>> new_pkgs;
            for (size_t i = 0; i < pkgs.size(); i+= 2) {
                vector<T> merged_sym(pkgs[i].first);
                merged_sym.insert(merged_sym.end(), pkgs[i+1].first.begin(), pkgs[i+1].first.end());

                u32 merge_val = pkgs[i].second + pkgs[i+1].second;
                new_pkgs.push_back(make_pair(merged_sym, merge_val));
            }
            pkgs = vector<pair<vector<T>, u32>>(original_freqs);
            pkgs.insert(pkgs.begin(), new_pkgs.begin(), new_pkgs.end());
            sort(pkgs.begin(), pkgs.end(), [](const pair<vector<T>, u32>& a, const pair<vector<T>, u32>& b) -> bool {
                return a.second < b.second;
            });
        }

        unordered_map<T, u32> map;
        for (size_t i = 0; i < num_pick_items; i++) {
            auto v = pkgs[i].first;
            for (T sym: v) {
                if (map.find(sym) == map.end()) {
                    map[sym] = 1;
                } else {
                    map[sym]++;
                }
            }
        }

        return map;
    }
    static vector<u32> generate_from_freqs(vector<u32> freqs) {
        vector<u32> res(freqs.size());

        //build priority queue
        priority_queue<Node*, vector<Node*>, Compare> q;
        for (u32 i = 0; i < freqs.size(); i++) {
            if (freqs[i] > 0) {
                Node* node = new Node(i, freqs[i]);
                node->left = nullptr;
                node->right = nullptr;
                q.push(node);
            }
        }

        Node *root = nullptr;

        if (q.size() == 1) {
            auto node = q.top();
            res[node->symbol] = 1;
            return res;
        }

        while (q.size() > 1) {
            Node *x = q.top();
            q.pop();
            Node * y = q.top();
            q.pop();

            Node* newNode = new Node(-1, x->value + y->value);
            newNode->left = x;
            newNode->right = y;
            root = newNode;
            q.push(newNode);
        }

        tree_traverse(root, res);

        return res;
    }

    static void tree_traverse(Node *root, vector<u32>& result) {
        do_tree_traverse(root, result, 0);
    }

    static void do_tree_traverse(Node *root, vector<u32>& result, u32 code_length) {
        if (root == nullptr) {
            return;
        }
        if (root->left == nullptr && root->right == nullptr) {
            result[root->symbol] = code_length;
            return;
        }
        do_tree_traverse(root->left, result, code_length + 1);
        do_tree_traverse(root->right, result, code_length + 1);
    }
};

//A simple downscaling algorithm using averaging.
inline std::vector<std::vector<unsigned char> > scale_down(std::vector<std::vector<unsigned char> > source_image, unsigned int source_width, unsigned int source_height, int factor){

    unsigned int scaled_height = (source_height+factor-1)/factor;
    unsigned int scaled_width = (source_width+factor-1)/factor;

    //Note that create_2d_vector automatically initializes the array to all-zero
    auto sums = create_2d_vector<unsigned int>(scaled_height,scaled_width);
    auto counts = create_2d_vector<unsigned int>(scaled_height,scaled_width);

    for(unsigned int y = 0; y < source_height; y++)
        for (unsigned int x = 0; x < source_width; x++){
            sums.at(y/factor).at(x/factor) += source_image.at(y).at(x);
            counts.at(y/factor).at(x/factor)++;
        }

    auto result = create_2d_vector<unsigned char>(scaled_height,scaled_width);
    for(unsigned int y = 0; y < scaled_height; y++)
        for (unsigned int x = 0; x < scaled_width; x++)
            result.at(y).at(x) = (unsigned char)((sums.at(y).at(x)+0.5)/counts.at(y).at(x));
    return result;
}

// Encode huffman a 8x8 quantized block
inline void output_block(OutputBitStream &output_stream, const vector<vector<int>>& block) {
    unordered_map<int, int> freqs_map;
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int e = block.at(y).at(x);
            if (freqs_map.find(e) == freqs_map.end()) {
                freqs_map[e] = 1;
            } else {
                freqs_map[e]++;
            }
        }
    // To vector of freqs


    //
    auto lengths_map = HuffmanCodingBuilder().package_merge(freqs_map, 16);
    int est_size = 0;
    int max_len = 0;
    for (auto l: lengths_map) {
        est_size += freqs_map[l.first] * l.second;
        if (l.second > max_len) {
            max_len = l.second;
        }
    }
    est_size /= 8;
    unordered_map<int, pair<u16, u16>> code_table = code_lengths_to_code_table(lengths_map);
    //Since block is 8x8, we are sure num of symbols less than 64 chars, wherese 1 byte = 256 chars.
    output_stream.push_byte((unsigned char)code_table.size());
    //
    for (auto c: code_table) {
        //Symbol
        output_stream.push_bits(c.first, 32);
        // Code length
        output_stream.push_byte(c.second.first);
    }

    int block_encoded_size = 0;
    int z_count = 0;
    for (int y = 0; y < 8; y++)
        for (int x = 0; x < 8; x++) {
            int e = block.at(y).at(x);
            if (e == 0) {
                z_count++;
            }
            auto code_pair = code_table[e];
            output_stream.push_bits(reverse_bits(code_pair.second, code_pair.first), code_pair.first);
            block_encoded_size += code_pair.first;
    }
}

template <typename T>
inline void output_plane(OutputBitStream &output_stream,
                         const vector<vector<int>>& quantize_vector,
                         vector<vector<T>>& in,
                         const vector<vector<Vector2D>> motion_vectors = vector<vector<Vector2D>>(),
                         bool is_P_frame = false,
                         int plane_scale = 1) {
    int in_h = (int)in.size();
    int in_w = (int)in[0].size();
    // each block is a vector of pairs of run size and its value
    vector<vector<pair<RunSize, int>>> data;
    for (int bh = 0; bh < in_h / 8; bh++)
        for (int bw = 0; bw < in_w / 8; bw++) {
            //What's this macro block index?
            if (is_P_frame) {
                int mb_x = bw * plane_scale / 2;
                int mb_y = bh * plane_scale / 2;
                if (mb_y < motion_vectors.size() && mb_x < motion_vectors[0].size()) {
                    Vector2D mv = motion_vectors.at(mb_y).at(mb_x);
                    if (mv.skip) {
                        //skip this block
                        continue;
                    }
                } else {
                    continue;
                }
            }

            //grab the 8x8 block starting at top left bh, bw
            auto block = create_2d_vector<double>(8, 8);
            for (unsigned int y = 0; y < 8; y++)
                for (unsigned int x = 0; x < 8; x++) {
                    block.at(y).at(x) = (double)in.at(bh * 8 + y).at(bw * 8 + x);
                }
            dct8_2d_transform(block);
            auto q_output = create_2d_vector<int>(8, 8);
            quantize(block, quantize_vector, q_output);

            //internal decompress to compute P-frames based on this decompressed frame
            {
                auto de_q_output = create_2d_vector<double>(8, 8);
                dequantize(q_output, quantize_vector, de_q_output);
                dct8_2d_inverse_transform(de_q_output);
                for (unsigned int y = 0; y < 8; y++)
                    for (unsigned int x = 0; x < 8; x++) {
                        in.at(bh * 8 + y).at(bw * 8 + x) = de_q_output.at(y).at(x);
                }
            }

            vector<int> flatten = zigzag_flatten(q_output);
            //Checking
            auto inflate = create_2d_vector<int>(8, 8);
            zigzag_inflate(flatten, inflate);
            assert(inflate == q_output);
            //End checking
            vector<pair<RunSize, int>> rle = rle_encoding(flatten);

            data.push_back(rle);
        }

    // Count freqs
    unordered_map<RunSize, int> freqs_map;
    for (unsigned int b = 0; b < data.size(); b++) {
        for (auto rs_pair: data[b]) {
            if (freqs_map.find(rs_pair.first) != freqs_map.end()) {
                freqs_map[rs_pair.first]++;
            } else {
                freqs_map[rs_pair.first] = 1;
            }
        }
    }

    auto lengths_map = HuffmanCodingBuilder().package_merge(freqs_map, 16);
    auto code_table = code_lengths_to_code_table(lengths_map);

    output_stream.push_byte((unsigned char)code_table.size());
    //
    for (auto c: code_table) {
        //Symbol
        u8 run = (u8)c.first.run;
        u8 size = (u8)c.first.size;
        output_stream.push_byte(run);
        output_stream.push_byte(size);
        // Code length
        output_stream.push_byte(c.second.first);
    }

    //Encode out data
    for (unsigned int b = 0; b < data.size(); b++) {
        for (auto rs_pair: data[b]) {
            auto code_pair = code_table[rs_pair.first];
            //output code
            output_stream.push_bits(reverse_bits(code_pair.second, code_pair.first), code_pair.first);
            //output value
            //e.g. min of 3 bits is -7
            // min of 5 bits is -31
            int min_val_in_bits_length = -((1<<rs_pair.first.size) - 1);
            if (rs_pair.second > 0) {
                output_stream.push_bits(rs_pair.second, rs_pair.first.size);
            } else {
                // if num is negative, we output the diff
                int diff = rs_pair.second - min_val_in_bits_length;
                output_stream.push_bits(diff, rs_pair.first.size);
            }
        }
    }
}

template <typename T>
inline void output(OutputBitStream& output_stream,
            const vector<vector<int>>& quantize_vector,
            vector<vector<T>>&Y,
            vector<vector<T>>&Cb,
            vector<vector<T>>&Cr) {
    output_plane(output_stream, quantize_vector, Y);
    output_plane(output_stream, quantize_vector, Cb);
    output_plane(output_stream, quantize_vector, Cr);
}

template <typename T>
inline void output_P_frame(OutputBitStream& output_stream,
                           const vector<vector<int>>& quantize_vector,
                           vector<vector<T>>&Y,
                           vector<vector<T>>&Cb,
                           vector<vector<T>>&Cr,
                           const vector<vector<Vector2D>>& motion_vects) {
    output_plane(output_stream, quantize_vector, Y, motion_vects, true, 1);
    output_plane(output_stream, quantize_vector, Cb, motion_vects, true, 2);
    output_plane(output_stream, quantize_vector, Cr, motion_vects, true, 2);
}


inline void compress_frame(YUVFrame420 &frame, YUVFrame420 &prev_frame, u32 width, u32 height, OutputBitStream &output_stream, const vector<vector<int>>& quantize_vector, bool is_intra) {
    u32 padded_height = next_mul(height, 8);
    u32 padded_width = next_mul(width, 8);

    auto Y = create_2d_vector<unsigned char>(padded_height, padded_width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Y.at(y).at(x) = frame.Y(x, y);

    //Extract the Cb plane into its own array
    int padded_c_h = next_mul((height + 1) / 2, 8);
    int padded_c_w = next_mul((width + 1) / 2, 8);
    auto Cb = create_2d_vector<unsigned char>(padded_c_h,padded_c_w);
    for(unsigned int y = 0; y < height/2; y++)
        for (unsigned int x = 0; x < width/2; x++)
            Cb.at(y).at(x) = frame.Cb(x, y);

    //Extract the Cr plane into its own array
    auto Cr = create_2d_vector<unsigned char>(padded_c_h,padded_c_w);
    for(unsigned int y = 0; y < height/2; y++)
        for (unsigned int x = 0; x < width/2; x++)
            Cr.at(y).at(x) = frame.Cr(x, y);

    //Output the 3 planes
    output(output_stream, quantize_vector, Y, Cb, Cr);
    output_stream.flush_to_byte();

    frame.set_Y_data(Y);
    frame.set_Cb_data(Cb);
    frame.set_Cr_data(Cr);
}


inline void compress_P_frames(YUVFrame420 &frame, YUVFrame420 &prev_frame, u32 width, u32 height, OutputBitStream &output_stream, const vector<vector<int>>& quantize_vector) {
    u32 padded_height = next_mul(height, 8);
    u32 padded_width = next_mul(width, 8);
    auto motion_vects = create_2d_vector<Vector2D>(height / 16, width / 16);
    //motion compensation

    auto Y = create_2d_vector<double>(padded_height, padded_width);

    //Extract the Cb plane into its own array
    int padded_c_h = next_mul((height + 1) / 2, 8);
    int padded_c_w = next_mul((width + 1) / 2, 8);
    //Extract the Cr plane into its own array
    auto Cb = create_2d_vector<double>(padded_c_h,padded_c_w);
    auto Cr = create_2d_vector<double>(padded_c_h,padded_c_w);

    for (unsigned int i = 0; i < height / 16; i++) {
        for (unsigned int j = 0; j < width / 16; j++) {
            Vector2D mv = block_match(frame, prev_frame, j * 16, i * 16);
            motion_vects.at(i).at(j) = mv;
            for (int m = 0; m < 16; m++)
                for (int n = 0; n < 16; n++) {
                    int x = j * 16 + n;
                    int y = i * 16 + m;
                    Y.at(y).at(x) = frame.Y(x, y) - prev_frame.Y(mv.x, mv.y);
                    Cb.at(y/2).at(x/2) = frame.Cb(x/2, y/2) - prev_frame.Cb(mv.x/2, mv.y/2);
                    Cr.at(y/2).at(x/2) = frame.Cr(x/2, y/2) - prev_frame.Cr(mv.x/2, mv.y/2);
//                    frame.Cb(x, y) =
//                    frame.Cb(x, y) -= prev_frame.Cb(mv.x, mv.y);
//                    frame.Cr(x, y) -= prev_frame.Cr(mv.x, mv.y);
                }
        }
    }
//
    for (unsigned int i = 0; i < height; i++) {
        for (unsigned int j = 0; j < width; j++) {
            if (i >= height / 16 * 16 || j >= width / 16 * 16) {
                Y.at(i).at(j) = frame.Y(j, i);
                Cb.at(i/2).at(j/2) = frame.Cb(j/2, i/2);
                Cr.at(i/2).at(j/2) = frame.Cr(j/2, i/2);
            }
        }
    }

    //Real difff
//    cerr<<"Real diff"<<endl;
//    for (int i = 0; i < 8; i++) {
//        for (int j = 0; j < 8; j++) {
//            cerr<<frame.Y(j, i)-prev_frame.Y(j, i)<<" ";
//        }
//        cerr<<endl;
//    }

    // Encode motion vectors for this frame
    // The number of vectors can always be calculated by the decoder
    // = (width/16) * (height/16) = number of macro blocks
    for (unsigned int y = 0; y < motion_vects.size(); y++) {
        for (unsigned int x = 0; x < motion_vects[0].size(); x++) {
            Vector2D mv = motion_vects.at(y).at(x);
            if (mv.skip) {
                output_stream.push_bit(0);
            } else {
                output_stream.push_bit(1);
                output_stream.push_u16((u16)mv.x);
                output_stream.push_u16((u16)mv.y);
            }
        }
    }
    // actual compress
    //prepare differences

    //Output the 3 planes
    output_P_frame(output_stream, quantize_vector, Y, Cb, Cr, motion_vects);
    output_stream.flush_to_byte();

//    cerr<<endl;
//    for (int i = 0; i < 8; i++) {
//        for (int j = 0; j < 8; j++) {
//            cerr<<Y.at(i).at(j)<<" ";
//        }
//        cerr<<endl;
//    }

    frame.set_Y_data(Y, prev_frame);
    frame.set_Cb_data(Cb, prev_frame);
    frame.set_Cr_data(Cr, prev_frame);

}

void do_compress(const char *argv[]) {
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
}

#endif
