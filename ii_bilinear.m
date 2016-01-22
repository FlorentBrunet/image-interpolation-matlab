function vals = ii_bilinear(img, x, y, nthreads, default_value)
%II_BILINEAR Bilinear interpolation.
%
% vals = ii_bilinear(img, x, y)
% vals = ii_bilinear(img, x, y, nthreads)
% vals = ii_bilinear(img, x, y, [], default_value)
% vals = ii_bilinear(img, x, y, nthreads, default_value)
%
% DESCRIPTION
%  "ii_bilinear" is a fast implementation of the standard bilinear
%  interpolation of an image. Parallelism is implemented using OpenMP. It
%  handles images with multiple channels. It also checks for points outside
%  of the image domain.
%
% INPUT ARGUMENTS
%  - img [array]: array of double values with a maximum of 3 dimensions.
%
%  - x, y [arrays]: locations where the image must be interpolated. They
%     can have any dimension. "x" and "y" must have the same number of 
%     elements. For best performances, try to use values in "x" and "y" 
%     that are consistent with the storage of the image.
%     Note: "x" are the abscissas (i.e. it corresponds to the columns of 
%     the image) and "y" are the ordinates (i.e. the rows). The upper left
%     corner of the image definition domain is (1,1).
%
%  - default_value [scalar, optional, default=0]: interpolated value when
%     the location is out of the image domain.
%
%  - nthreads [scalar, optional]: default number of threads used to compute
%     the interpolation. This default value is defined by the C++ macro
%     DEFAULT_NTHREADS in "ii_bilinear.cpp".
%
% OUTPUT ARGUMENTS
%  - vals [array]: interpolated values. If "dims" are the dimensions of the
%     input "x" then the dimension of "vals" is [dims nchannels] where
%     "nchannels" is the number of color channel of the input image "img".
%
% 2012/07/17, Florent Brunet

% This file is part of II.
%
% II is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% II is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with II.  If not, see <http://www.gnu.org/licenses/>.
