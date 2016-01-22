function test_compare
%TEST_COMPARE Simple demo/test comparing the ii_bicubic and the ii_bilinear functions.
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

close all;

% Load an image
img = double(imread('lena.jpg'))/255;

% Interpolate the image using "ii_bicubic"
[x y] = meshgrid(1:0.5:512);
interp_img_linear = ii_bilinear(img, x, y);
interp_img_cubic = ii_bicubic(img, x, y);

figure(1);
imshow(img);
title('Original');

figure(2);
imshow(interp_img_linear);
title('Bilinear interpolation');

figure(3);
imshow(interp_img_cubic);
title('Bicubic interpolation');

figure(4);
imagesc(mean(interp_img_cubic - interp_img_linear, 3));
daspect([1 1 1]);
colorbar;
title('Difference between bilinear and bicubic');
