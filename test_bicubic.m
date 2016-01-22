function test_bicubic
%TEST_BICUBIC Simple demo/test of the ii_bicubic function.
%
% 2011/05/23, Florent Brunet

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

close all;

% Load an image
img = double(imread('lena.jpg'))/255;

% Interpolate the image using "ii_bicubic"
[x y] = meshgrid(1:0.5:512);
interp_img = ii_bicubic(img, x, y);

figure(1);
imshow(img);

figure(2);
imshow(interp_img);
