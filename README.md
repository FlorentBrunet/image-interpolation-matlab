# image-interpolation-matlab

## Description

image-interpolation-matlab is a small Matlab toolbox that provides really fast and convenient image interpolation routines. It is fast since it relies on native binary code (Mex-files) and it is implemented in parallel. It is convenient because it can handle images with multiple channels (contrarily to the "interp2" function of Matlab). It supports bicubic and bilinear interpolation schemes. I hope one day to have some time to implement other algorithms.

## Instructions

If you download the source code of this toolbox (which is the preferred way), you will need to compile the mex-files. To do so, just run the script named "ii_compile_and_setup.m". There are a few options with self-explaining names at the beginning of this script if you want to tune the toolbox according to your needs and/or your hardware.

