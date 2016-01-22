/* ii_bilinear.cpp
 * 2012/07/17, Florent Brunet
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

// This function handles special cases (borders, out of domain, etc.)
double other_cases(const double *img, int height, int width, const double& default_value,
                   int rx, int ry,
                   const double& cx, const double& cy) {

    if ((rx == width-1) && (ry < height-1)) {
        // Right border
        
        if (cx > 0) {
            return default_value;
        } else {
            double v00 = img[rx * height + ry];
            double v01 = img[rx * height + ry+1];
            return v00 * (1.0 - cy) + v01 * cy;
        }
        
    } else if ((rx < width-1) && (ry == height-1)) {
        // Bottom border
        
        if (cy > 0) {
            return default_value;
        } else {
            double v00 = img[rx * height + ry];
            double v10 = img[(rx+1) * height + ry];
            return v00 * (1.0 - cx) + v10 * cx;
        }
        
    } else if ((rx == width-1) && (ry == height-1)) {
        // Bottom right pixel
        
        if ((cx > 0) || (cy > 0)) {
            return default_value;
        } else {
            return img[rx * height + ry];
        }
        
    } else {
        // Default value out of the bounds
        return default_value;
    }
}


void bilinear(const double *img, int height, int width, int nchannels,
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
            double omcx = 1.0 - cx;
            rx--;

            double cy = y[i];
            int ry = static_cast<int>(cy);
            cy -= ry;
            ry--;

            if ((rx >= 0) && (rx < width-1) && (ry >= 0) && (ry < height-1)) {
                // This is the general case (i.e. the 'middle' of the image)
                // The other cases are deported in the function 'other_cases'
                // (it makes the inner code of the loop smaller and consequently
                // reduces the number of 'instruction cache misses')
                const double *cur_chan = &(chan[rx*height + ry]);
                
                double v00 = *cur_chan;
                cur_chan++;
                double v01 = *cur_chan;
                cur_chan += height;
                double v11 = *cur_chan;
                cur_chan--;
                double v10 = *cur_chan++;

                interp_chan[i] = (v00 * omcx + v10 * cx) * (1.0 - cy) + (v01 * omcx + v11 * cx) * cy;
            } else
                interp_chan[i] = other_cases(chan, height, width, default_value, rx, ry, cx, cy);
        }
    }
}

// ii_bilinear(img, x, y, [nthreads], [default])
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
    
    bilinear(img, (int)height, (int)width, (int)nchannels, x, y, (int)npts, interp, nthreads, default_value);
}
