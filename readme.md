This program is a simple implementation of MPEG-1 video compression standard, which implements I-frame and P-frame.

# Sample/Test video files
Sample video files can be downloaded from [Xiph](https://media.xiph.org/video/derf/).

# Usage
- This implementation only works with raw yuv420 files, which means that you'll need `ffmpeg` installed to convert between standard video files to raw yuv420 files.
- To compress a raw yuv420 file, convert a standard video file to raw yuv420 first:

`ffmpeg -i source_video.y4m -f rawvideo -pixel_format yuv420p - > input_video.raw`,

then compress the raw yuv420 file:

`./pmpeg 352 288 medium < input_video.raw > compressed.vid`

where 352 and 288 are width and height of the source video, medium is compress quality, quality options are low/medium/high.

- To decompress a compressed video file:

`./pmpeg -d < compressed.vid > output_video.raw`,

then convert the raw yuv420 file to y4m by 

`ffmpeg -f rawvideo -pixel_format yuv420p -framerate 30 -video_size 352x288 -i - -f yuv4mpegpipe output_video.y4m < output_video.raw`

remember to input correct video dimensions.

- You can play the output y4m video file with `ffplay output_video.y4m` to see the quality of the compressed video.

# Compression scheme
- An I-frame is technically compressed with following transformations and encodings:
    - DCT
    - Quantization
    - RLE on the zigzag re-ordered coefficients
    - Dynamic Huffman coding of run-size symbols (I did not differentiate AC and DC coefficients).
- With a P-frame, motion vectors are computed first, then the differences between the current frame and the previous decompressed frame are calculated. Then the same set of transformations which are used in I-frames is performed.

# Bitstream format
The bitstream format of a compressed file contains a header and compressed frames. A header is as follows:
- Video width: 4 bytes
- Video height: 4 bytes
- Quality: 2 bits (0 = LOW, 1 = MEDIUM, 2 = HIGH)
Then there's 1 byte flag indicates whether there is a next frame. The decoder stops when reaching a 0 flag. Followed by this flag is a compressed frame, either I-frame or P-frame.
An I-frame consists of:
- The number of Huffman symbols N: 1 byte
- N code words, each code word takes 4 bytes for symbol value, 1 byte symbol length, so 5 bytes per symbol.
- Run-size symbols of compressed blocks:
    - Code value
    - Value of the cofficient (if the coef != zero)
    (for example, if a run size (3, 3)'s code is 10001, the coef is 7 (3 bits) then we output 10001 111,
    but if the coef is -7 then we output 10001 000)
A P-frame consists the same bitstream as an I-frame, but before that we output `width / 16 * height / 16` number of motion vectors. More specfically,
- Has Motion flag: 1 bit, if this flag is zero, the decoder skips this block, if this block has a motion vector, we encode next the motion vector value
- Motion vector's x: 2 bytes
- Motion vector's y: 2 bytes

# Implementation notes
- The output video frame sequence (raw YUV420) contains groups of pictures (GOPs), where each GOP contains 4 frames: 1 I-frame and 3 P-frames. The first frame of the output video frame sequence is an I-frame.
- Each frame is divided into 8x8 blocks, each block then is passed through DCT, zigzag-reordering quantization, and then a dynamic Huffman coding (similar to JPEG compression). A fast version of DCT is adopted from this work of Nayuki (https://www.nayuki.io/page/fast-discrete-cosine-transform-algorithms). This fast DCT performs 6-8 time faster than a naive DCT. The output matrix of 2d DCT is then quantized, flattened into a 1D array and then sorted using a zigzag matrix as described in JPEG specs. This zigzag-reordering is expected to group zero coefficients to the end of the 1D array. Run Length Encoding then encodes the 1D array coefs into run-size symbols. Each symbol is a struct RunSize, which indicates a run (the number of zeroes before this symbol), and a size (the number of bits needed to encode the coef). By only concerning about coef size (number of bits), we can reduce the number of symbols used by Huffman coding (otherwise there would be thousands of symbols).
- For I-frames, the whole frame is encoded so that the decoder doesn't need further information to decode the frame.
- For P-frames, we encode the differences between the current frame and the previous decompressed frame (it can be an I-frame or a P-frame) to reduce prediction errors from compression (actually quantization). For each block, we implement a not-so-naive block matching algorithm called Three Step Search to find motion vectors that can best match the current block. Mean absolute difference (MAD) metric is used to compare the errors between the blocks. To speed up this motion compensation step, we compute motion vectors for 16x16 macro blocks(a 16x16 macro block contains 4 Y blocks, 1 Cb block and 1 Cr block). If the motion vector found is (0, 0) and the MAD is below a THRESHOLD of 3.0, we SKIP the block and notify the decoder to use the block at the same position from the previous frame. This SKIP flag only needs 1 bit so that it significantly reduces the compressed video size. We encode motion vector components (x, y) as 16 bit unsigned integers (so the program only supports each video dimension upto 65536 pixels).

# Compresion ratio
Some compression ratios of sample videos are as follows (raw y4m files are ~45.6MB)
- Harbour:
    - high quality: ~20 (2.252MB)
- News:
    - high quality: ~31.82 (1.442 MB)
- Flower:
    - high quality: ~18.58 (2.454 MB)
