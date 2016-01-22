/* ii_bicubic.cpp
 * 2011/05/23, Florent Brunet
 *
 * This file is part of II.
 *
 * II is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * II is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with II.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mex.h"
#include <omp.h>

#ifndef DEFAULT_NTHREADS
#define DEFAULT_NTHREADS 8
#endif

// This function handles all the special cases (borders, out of domain, etc.)
double other_cases(const double *img, int height, int width, const double& default_value,
                 int rx, int ry,
                 const double& cx, const double& cx2, const double& cx3,
                 const double& cy, const double& cy2, const double& cy3) {
    double p0 = 2*cy2-cy3-cy;
    double p1 = 3*cy3-5*cy2+2;
    double p2 = 4*cy2-3*cy3+cy;
    double p3 = cy3-cy2;
    
    if ((rx == 0) && (ry == 0)) {
        // CASE 1 ---------------------------------------------------------

        double v10 = 3*img[rx*height + ry] - 3*img[rx*height + ry+1] + img[rx*height + ry+2];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        double v20 = 3*img[(rx+1)*height + ry] - 3*img[(rx+1)*height + ry+1] + img[(rx+1)*height + ry+2];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry+2];

        double v30 = 3*img[(rx+2)*height + ry] - 3*img[(rx+2)*height + ry+1] + img[(rx+2)*height + ry+2];
        double v31 = img[(rx+2)*height + ry];
        double v32 = img[(rx+2)*height + ry+1];
        double v33 = img[(rx+2)*height + ry+2];

        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;
        double v0 = 3*v1 - 3*v2 + v3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx == 0) && (ry >= 1) && (ry < height-2)) {
        // CASE 5 ---------------------------------------------------------

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        double v20 = img[(rx+1)*height + ry-1];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry+2];

        double v30 = img[(rx+2)*height + ry-1];
        double v31 = img[(rx+2)*height + ry];
        double v32 = img[(rx+2)*height + ry+1];
        double v33 = img[(rx+2)*height + ry+2];

        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;
        double v0 = 3*v1 - 3*v2 + v3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx == 0) && (ry == height-2)) {
        // CASE 9 ---------------------------------------------------------

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry-1] - 3*img[rx*height + ry] + 3*img[rx*height + ry+1];

        double v20 = img[(rx+1)*height + ry-1];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry-1] - 3*img[(rx+1)*height + ry] + 3*img[(rx+1)*height + ry+1];
        
        double v30 = img[(rx+2)*height + ry-1];
        double v31 = img[(rx+2)*height + ry];
        double v32 = img[(rx+2)*height + ry+1];
        double v33 = img[(rx+2)*height + ry-1] - 3*img[(rx+2)*height + ry] + 3*img[(rx+2)*height + ry+1];

        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;
        double v0 = 3*v1 - 3*v2 + v3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx == 0) && (ry == height-1) && (cy == 0)) {
        // CASE 13 --------------------------------------------------------

        double v1 = img[rx*height+ry];
        double v2 = img[(rx+1)*height+ry];
        double v3 = img[(rx+2)*height+ry];
        double v0 = 3*v1 - 3*v2 + v3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 2;
        
    } else if ((rx >= 1) && (rx < width-2) && (ry == 0)) {
        // CASE 2 ---------------------------------------------------------

        double v00 = 3*img[(rx-1)*height + ry] - 3*img[(rx-1)*height + ry+1] + img[(rx-1)*height + ry+2];
        double v01 = img[(rx-1)*height + ry];
        double v02 = img[(rx-1)*height + ry+1];
        double v03 = img[(rx-1)*height + ry+2];

        double v10 = 3*img[rx*height + ry] - 3*img[rx*height + ry+1] + img[rx*height + ry+2];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        double v20 = 3*img[(rx+1)*height + ry] - 3*img[(rx+1)*height + ry+1] + img[(rx+1)*height + ry+2];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry+2];

        double v30 = 3*img[(rx+2)*height + ry] - 3*img[(rx+2)*height + ry+1] + img[(rx+2)*height + ry+2];
        double v31 = img[(rx+2)*height + ry];
        double v32 = img[(rx+2)*height + ry+1];
        double v33 = img[(rx+2)*height + ry+2];

        double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx >= 1) && (rx < width-2) && (ry == height-2)) {
        // CASE 10 --------------------------------------------------------
        
        double v00 = img[(rx-1)*height + ry-1];
        double v01 = img[(rx-1)*height + ry];
        double v02 = img[(rx-1)*height + ry+1];
        double v03 = img[(rx-1)*height + ry-1] - 3*img[(rx-1)*height + ry] + 3*img[(rx-1)*height + ry+1];
        
        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry-1] - 3*img[rx*height + ry] + 3*img[rx*height + ry+1];

        double v20 = img[(rx+1)*height + ry-1];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry-1] - 3*img[(rx+1)*height + ry] + 3*img[(rx+1)*height + ry+1];
        
        double v30 = img[(rx+2)*height + ry-1];
        double v31 = img[(rx+2)*height + ry];
        double v32 = img[(rx+2)*height + ry+1];
        double v33 = img[(rx+2)*height + ry-1] - 3*img[(rx+2)*height + ry] + 3*img[(rx+2)*height + ry+1];

        double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx >= 1) && (rx < width-2) && (ry == height-1) && (cy == 0)) {
        // CASE 14 --------------------------------------------------------
        
        double v0 = img[(rx-1)*height + ry];
        double v1 = img[rx*height + ry];
        double v2 = img[(rx+1)*height + ry];        
        double v3 = img[(rx+2)*height + ry];

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 2;
        
    } else if ((rx == width-2) && (ry == 0)) {
        // CASE 3 ---------------------------------------------------------

        double v00 = 3*img[(rx-1)*height + ry] - 3*img[(rx-1)*height + ry+1] + img[(rx-1)*height + ry+2];
        double v01 = img[(rx-1)*height + ry];
        double v02 = img[(rx-1)*height + ry+1];
        double v03 = img[(rx-1)*height + ry+2];

        double v10 = 3*img[rx*height + ry] - 3*img[rx*height + ry+1] + img[rx*height + ry+2];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        double v20 = 3*img[(rx+1)*height + ry] - 3*img[(rx+1)*height + ry+1] + img[(rx+1)*height + ry+2];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry+2];

        double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = 3*v2 - 3*v1 + v0;
        
        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
    
    } else if ((rx == width-2) && (ry >= 1) && (ry < height-2)) {
        // CASE 7 ---------------------------------------------------------

        double v00 = img[(rx-1)*height + ry-1];
        double v01 = img[(rx-1)*height + ry];
        double v02 = img[(rx-1)*height + ry+1];
        double v03 = img[(rx-1)*height + ry+2];

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        double v20 = img[(rx+1)*height + ry-1];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry+2];

        double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = 3*v2 - 3*v1 + v0;
        
        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx == width-2) && (ry == height-2)) {
        // CASE 11 --------------------------------------------------------

        double v00 = img[(rx-1)*height + ry-1];
        double v01 = img[(rx-1)*height + ry];
        double v02 = img[(rx-1)*height + ry+1];
        double v03 = img[(rx-1)*height + ry-1] - 3*img[(rx-1)*height + ry] + 3*img[(rx-1)*height + ry+1];

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry-1] - 3*img[rx*height + ry] + 3*img[rx*height + ry+1];

        double v20 = img[(rx+1)*height + ry-1];
        double v21 = img[(rx+1)*height + ry];
        double v22 = img[(rx+1)*height + ry+1];
        double v23 = img[(rx+1)*height + ry-1] - 3*img[(rx+1)*height + ry] + 3*img[(rx+1)*height + ry+1];

        double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
        double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
        double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
        double v3 = 3*v2 - 3*v1 + v0;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
        
    } else if ((rx == width-2) && (ry == height-1) && (cy == 0)) {
        // CASE 15 --------------------------------------------------------

        double v0 = img[(rx-1)*height + ry];
        double v1 = img[rx*height + ry];
        double v2 = img[(rx+1)*height + ry];
        double v3 = 3*v2 - 3*v1 + v0;

        return (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 2;
        
    } else if ((rx == width-1) && (cx == 0) && (ry == 0)) {
        // CASE 4 ---------------------------------------------------------

        double v10 = 3*img[rx*height + ry] - 3*img[rx*height + ry+1] + img[rx*height + ry+2];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        return (v10*p0 + v11*p1 + v12*p2 + v13*p3) / 2;
        
    } else if ((rx == width-1) && (cx == 0) && (ry >= 1) && (ry < height-2)) {
        // CASE 8 ---------------------------------------------------------

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = img[rx*height + ry+2];

        return (v10*p0 + v11*p1 + v12*p2 + v13*p3) / 2;
        
    } else if ((rx == width-1) && (cx == 0) && (ry == height-2)) {
        // CASE 12 --------------------------------------------------------

        double v10 = img[rx*height + ry-1];
        double v11 = img[rx*height + ry];
        double v12 = img[rx*height + ry+1];
        double v13 = 3*v12 - 3*v11 + v10;

        return (v10*p0 + v11*p1 + v12*p2 + v13*p3) / 2;
        
    } else if ((rx == width-1) && (cx == 0) && (ry == height-1) && (cy == 0)) {
        // CASE 16 --------------------------------------------------------
        
        return img[rx*height+ry];
                
    } else
        // Default value out of the bounds
        return default_value;
}


// Actuel bicubic interpolation
void bicubic(const double *img, int height, int width, int nchannels,
             double *x, double *y, int npts,
             double *interp, int nthreads, const double& default_value) {

    for (int c = 0; c < nchannels; ++c) {
        const double *chan = img + c*height*width;
        double *interp_chan = interp + c*npts;
        
        #pragma omp parallel for num_threads(nthreads), schedule(guided)
        for (int i = 0; i < npts; ++i) {
            double cx = x[i];
            int rx = static_cast<int>(cx);
            cx -= rx;
            double cx2 = cx*cx;
            double cx3 = cx*cx2;
            rx--;

            double cy = y[i];
            int ry = static_cast<int>(cy);
            cy -= ry;
            double cy2 = cy*cy;
            double cy3 = cy*cy2;
            ry--;

            if ((rx >= 1) && (rx < width-2) && (ry >= 1) && (ry < height-2)) {
                // This is the general cases (i.e. the 'middle' of the image)
                // The other cases are deported in the function 'other_cases'
                // (it makes the inner code of the loop smaller and consequently
                // reduces the number of 'instruction cache misses')
                const double *cur_chan_base = &(chan[(rx-1)*height + ry-1]);
                const double *cur_chan = cur_chan_base;
                
                double v00 = *cur_chan++;
                double v01 = *cur_chan++;
                double v02 = *cur_chan++;
                double v03 = *cur_chan++;

                cur_chan_base += height;
                cur_chan = cur_chan_base;
                double v10 = *cur_chan++;
                double v11 = *cur_chan++;
                double v12 = *cur_chan++;
                double v13 = *cur_chan++;

                cur_chan_base += height;
                cur_chan = cur_chan_base;
                double v20 = *cur_chan++;
                double v21 = *cur_chan++;
                double v22 = *cur_chan++;
                double v23 = *cur_chan++;

                cur_chan = cur_chan_base + height;
                double v30 = *cur_chan++;
                double v31 = *cur_chan++;
                double v32 = *cur_chan++;
                double v33 = *cur_chan++;

                double p0 = 2*cy2-cy3-cy;
                double p1 = 3*cy3-5*cy2+2;
                double p2 = 4*cy2-3*cy3+cy;
                double p3 = cy3-cy2;

                double v0 = v00*p0 + v01*p1 + v02*p2 + v03*p3;
                double v1 = v10*p0 + v11*p1 + v12*p2 + v13*p3;
                double v2 = v20*p0 + v21*p1 + v22*p2 + v23*p3;
                double v3 = v30*p0 + v31*p1 + v32*p2 + v33*p3;

                interp_chan[i] = (v0*(2*cx2-cx3-cx) + v1*(3*cx3-5*cx2+2) + v2*(4*cx2-3*cx3+cx) + v3*(cx3-cx2)) / 4;
            } else
                interp_chan[i] = other_cases(chan, height, width, default_value, rx, ry, cx, cx2, cx3, cy, cy2, cy3);
        }
    }
}

// ii_bicubic(img, x, y, [nthreads], [default])
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if ((nrhs < 3) || (nrhs > 5))
        mexErrMsgTxt("The number of input arguments must be comprised between 3 and 5.");
    
    // prhs[0]: img
    if (!mxIsDouble(prhs[0]))
        mexErrMsgTxt("The type of the input image must be double.");
    
    double *img = mxGetPr(prhs[0]);
    mwSize ndims_img = mxGetNumberOfDimensions(prhs[0]);
    if (ndims_img > 3)
        mexErrMsgTxt("The number of dimensions of the input image cannot be greater than 3.");
    
    const mwSize *dims_img = mxGetDimensions(prhs[0]);
    size_t height = dims_img[0];
    size_t width = 1;
    if (ndims_img > 1)
        width = dims_img[1];
    size_t nchannels = 1;
    if (ndims_img > 2)
        nchannels = dims_img[2];
    if (height <= 3)
        mexErrMsgTxt("The height of the input image cannot be smaller than 3.");
    if (width <= 3)
        mexErrMsgTxt("The width of the input image cannot be smaller than 3.");
    
    // prhs[1]: x
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("The type of the input argument 'x' must be double.");
    double *x = mxGetPr(prhs[1]);
    size_t npts = mxGetNumberOfElements(prhs[1]);
    mwSize ndims_x = mxGetNumberOfDimensions(prhs[1]);
    const mwSize *dims_x = mxGetDimensions(prhs[1]);
    
    // prhs[2]: y
    if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("The type of the input argument 'y' must be double.");
    if (npts != mxGetNumberOfElements(prhs[2]))
        mexErrMsgTxt("'x' and 'y' must have the same number of elements.");
    double *y = mxGetPr(prhs[2]);
    
    // prhs[3]: nthreads
    int nthreads = DEFAULT_NTHREADS;
    if (nrhs > 3) {
        if (!mxIsEmpty(prhs[3])) {
            if (!mxIsNumeric(prhs[3]))
                mexErrMsgTxt("'nthreads' must be numeric.");
            nthreads = static_cast<int>(mxGetScalar(prhs[3]));
            if (nthreads < 1)
                mexWarnMsgTxt("'nthreads' should be greater than 1. Using the default number of threads.");
        }
    }
    
    // prhs[4]: default_value
    double default_value = 0.0;
    if (nrhs > 4) {
        if (!mxIsDouble(prhs[4]))
            mexErrMsgTxt("'default_value' must be of type double.");
        default_value = mxGetScalar(prhs[4]);
    }

    // plhs[0]: interpolated values
    mwSize ndims_out;
    mwSize *dims_out = 0;
    
    ndims_out = ndims_x;
    if (ndims_img == 3)
        ndims_out++;
    dims_out = new mwSize[ndims_out];
    for (int i = 0; i < ndims_x; ++i)
        dims_out[i] = dims_x[i];
    if (ndims_img == 3)
        dims_out[ndims_out-1] = dims_img[2];
    plhs[0] = mxCreateNumericArray(ndims_out, dims_out, mxDOUBLE_CLASS, mxREAL);
    double *interp = mxGetPr(plhs[0]);
    
    bicubic(img, (int)height, (int)width, (int)nchannels, x, y, (int)npts, interp, nthreads, default_value);
}
