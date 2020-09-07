//
//  pmpeg.cpp
//  A5
//
//  Created by TonyNguyen on 9/6/20.
//  Copyright © 2020 Phuc Nguyen. All rights reserved.
//

#include <iostream>
#include "compress.h"
#include "decompress.h"

using namespace std;

void print_usage(const char *argv[]) {
    cerr<<"Usage: "<<argv[0]<<" [-d] [low/medium/high] <input_file> <output_file>"<<endl;
    cerr<<"Compress: "<<argv[0]<<"width height low/medium/high < input_file > output_file"<<endl;
    cerr<<"Decompress: "<<argv[0]<<"-d input_file output_file"<<endl;
}

int main(int argc, const char *argv[]) {
    bool decompress_mode = false;
    for (unsigned int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0) {
            decompress_mode = true;
            break;
        }
    }

    if ((!decompress_mode || argc != 2) && argc != 4) {
        print_usage(argv);
        return -1;
    }

    if (decompress_mode) {
        do_decompress();
    } else {
        do_compress(argv);
    }

    return 0;
}
