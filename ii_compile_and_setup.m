function ii_compile_and_setup
%II_COMPILE_AND_SETUP Matlab script for compiling the II toolbox.
%
% IMPORTANT:
% Change the options at the beginning of this script in order to tune the
% compilation to your hardware and your system.
% Then launch this script from the folder of the II toolbox.
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


% Uncomment the line corresponding to your system:
system_name = 'linux';
% system_name = 'windows';

% Default number of threads. Usually, the actual number of core on your system
% is a good choice. You can experiment the toolbox with "test_timings" to
% determine the optimal number of threads.
% If you set "nthreads" to 0 then multithreading (and OpenMP) won't be used.
% Note that it is always possible to override this value when calling the
% functions of the toolbox.
nthreads = 8;
% nthreads = 0

delete *.obj
delete *.mex*
delete *.pdb
delete *.ilk

if strcmpi(system_name, 'linux')
    %% Linux
    if nthreads <= 0
        para_str = '-fopenmp -DDEFAULT_NTHREADS=1';
    else
        para_str = sprintf('-fopenmp -DDEFAULT_NTHREADS=%d', nthreads);
    end
        
    eval(sprintf('mex -cxx CXXFLAGS="\\$CXXFLAGS %s" LDFLAGS="\\$LDFLAGS -fopenmp" -v ii_bilinear.cpp', para_str));
    eval(sprintf('mex -cxx CXXFLAGS="\\$CXXFLAGS %s" LDFLAGS="\\$LDFLAGS -fopenmp" -v ii_bicubic.cpp', para_str));
    
elseif strcmpi(system_name, 'windows');
    %% Windows
    if nthreads <= 0
        para_str = '/openmp -DDEFAULT_NTHREADS=1';
    else
        para_str = sprintf('/openmp -DDEFAULT_NTHREADS=%d', nthreads);
    end
    
    eval(sprintf('mex -cxx COMPFLAGS="$COMPFLAGS %s" ii_bicubic.cpp', para_str));
    eval(sprintf('mex -cxx COMPFLAGS="$COMPFLAGS %s" ii_bilinear.cpp', para_str));
else
    error('Invalid system name');
end

% Add the toolbox to the Matlab path and save it
rmpath(pwd);
addpath(pwd);
savepath;
