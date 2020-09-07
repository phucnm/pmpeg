//
//  motion_compensation.hpp
//  A5
//
//  Created by TonyNguyen on 8/9/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#ifndef MOTION_COMPENSATION_HPP
#define MOTION_COMPENSATION_HPP

#include <iostream>
#include <vector>
#include "yuv_stream.hpp"

using namespace std;

const vector<int> x_offset{0, -1,  0,  1, -1, 1, -1, 0, 1};
const vector<int> y_offset{0, -1, -1, -1,  0, 0,  1, 1, 1};
const int MACRO_BLK_SIZE = 16;
const double ERROR_THRESHOLD = 3;

struct Vector2D {
    int x, y;
    bool skip = false;
};

// Compute MAD at the top left (x, y) with size of macro block size
inline double compute_AAD(YUVFrame420& frame, YUVFrame420& prev_frame, int block_x, int block_y, int ref_x, int ref_y) {
    double aad = 0;
    for (int i = 0; i < MACRO_BLK_SIZE; i++) {
        for (int j = 0; j < MACRO_BLK_SIZE; j++) {
            aad += abs(frame.Y(block_x + i, block_y + j) - prev_frame.Y(ref_x + i, ref_y + j));
        }
    }
    aad /= 256; // N^2 = 16*16
    return aad;
}

inline Vector2D do_3_step_search(YUVFrame420& frame, YUVFrame420 &prev_frame, int block_x, int block_y, int step_size) {
    Vector2D best_mv {block_x, block_y, .skip = true};
    double min_AAD = compute_AAD(frame, prev_frame, block_x, block_y, block_x, block_y);
    for (unsigned int i = 0; i < x_offset.size(); i++) {
        int mv_x = block_x + step_size * x_offset[i];
        int mv_y = block_y + step_size * y_offset[i];
        if (mv_x >=0 && mv_x < frame.get_width() - MACRO_BLK_SIZE &&
            mv_y >=0 && mv_y < frame.get_height() - MACRO_BLK_SIZE) {
            double AAD = compute_AAD(frame, prev_frame, block_y, block_y, mv_x , mv_y);
            if (AAD < min_AAD && AAD < ERROR_THRESHOLD) {
                min_AAD = AAD;
                best_mv = Vector2D{mv_x, mv_y};
            }
        }
    }

    // half way stop
    if (best_mv.x == block_x && best_mv.y == block_y) {
        best_mv.skip = true;
        return best_mv;
    }

    if (step_size == 1)
        return best_mv;

    return do_3_step_search(frame, prev_frame, best_mv.x, best_mv.y, step_size / 2);
}

inline Vector2D block_match(YUVFrame420& frame, YUVFrame420 &prev_frame, int block_x, int block_y) {
    return do_3_step_search(frame, prev_frame, block_x, block_y, 4);
}
#endif
