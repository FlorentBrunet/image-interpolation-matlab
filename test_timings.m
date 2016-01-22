function test_timings
%TEST_TIMINGS Several tests of timings.
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
clc;

comparison_with_matlab(4, 5, 512);
% optimal_number_of_threads

%% Compare II against Matlab
function comparison_with_matlab(nthreads, ntrials, sz)
% nthreads: number of threads for II
% ntrials: number of times to restart the experiment to get godd averaged results
% sz: size of the image

[x y] = meshgrid(linspace(1, sz, 2*sz-1));

t_lin_mat = 0;
t_cub_mat = 0;
t_lin_ii = 0;
t_cub_ii = 0;

err_lin = 0;
err_cub = 0;

fprintf('\n[II/Matlab] nthreads=%d, ntrials=%d, sz=%dx%d\n\n', nthreads, ntrials, sz, sz);

for i = 1:ntrials
    fprintf('[II/Matlab] Trial %d\n', i);
    
    img = rand(sz);
    
    % Matlab
    tic;
    val_img_lin_matlab = interp2(img, x, y, 'linear');
    t_lin_mat = t_lin_mat + toc;

    tic;
    val_img_cub_matlab = interp2(img, x, y, 'cubic');
    t_cub_mat = t_cub_mat + toc;


    % II
    tic;
    val_img_lin_ii = ii_bilinear(img, x, y, nthreads);
    t_lin_ii = t_lin_ii + toc;

    tic;
    val_img_cub_ii = ii_bicubic(img, x, y, nthreads);
    t_cub_ii = t_cub_ii + toc;
    
    % Errors
    err_lin = err_lin + norm(val_img_lin_matlab - val_img_lin_ii);
    err_cub = err_cub + norm(val_img_cub_matlab - val_img_cub_ii);
    
    clear('img');
end

fprintf('\n[II/Matlab] Error bilinear: %f\n', err_lin/ntrials);
fprintf('[II/Matlab] Error bicubic: %f\n', err_cub/ntrials);

fprintf('\n[II/Matlab] Timing bilinear Matlab: %f s\n', t_lin_mat/ntrials);
fprintf('[II/Matlab] Timing bilinear II    : %f s\n', t_lin_ii/ntrials);
fprintf('[II/Matlab] Bilinear: II %f times faster than Matlab\n', t_lin_mat/t_lin_ii);

fprintf('\n[II/Matlab] Timing bicubic Matlab: %f s\n', t_cub_mat/ntrials);
fprintf('[II/Matlab] Timing bicubic II    : %f s\n', t_cub_ii/ntrials);
fprintf('[II/Matlab] Bicubic: II %f times faster than Matlab\n', t_cub_mat/t_cub_ii);

fprintf('\n[II/Matlab] II_bilin %f times faster than II_bicub\n', t_cub_ii/t_lin_ii);


%% Show the efficiency of II in function of the number of threads
function optimal_number_of_threads
number_of_threads_to_test = 1:16;
ntrials = 500;
sz = 512;

[x y] = meshgrid(linspace(1, sz, 2*sz-1));

t_lin = zeros(ntrials, numel(number_of_threads_to_test));
t_cub = zeros(ntrials, numel(number_of_threads_to_test));

fprintf('[II # threads]\n\n');

for i_nthreads = 1:numel(number_of_threads_to_test)
    nthreads = number_of_threads_to_test(i_nthreads);
    
    for i_trial = 1:ntrials
        fprintf('[II # threads] nthreads=%d, trial=%d\n', nthreads, i_trial);
        
        img = rand(sz);
        
        tic;
        ii_bilinear(img, x, y, nthreads);
        t_lin(i_trial, i_nthreads) = toc;
        
        tic;
        ii_bicubic(img, x, y, nthreads);
        t_cub(i_trial, i_nthreads) = toc;
        
        clear('img');
    end
end


% "Clean" the results (i.e. throw away the 10% best and 10% worst timings)
neject = max(0, round(0.1*ntrials));

t_lin_clean = zeros(ntrials-2*neject, numel(number_of_threads_to_test));
t_cub_clean = zeros(ntrials-2*neject, numel(number_of_threads_to_test));
for i_nthreads = 1:numel(number_of_threads_to_test)
    col = t_lin(:,i_nthreads);
    col = sort(col);
    col(1:neject) = [];
    col(end-neject+1:end) = [];
    t_lin_clean(:,i_nthreads) = col;
    
    col = t_cub(:,i_nthreads);
    col = sort(col);
    col(1:neject) = [];
    col(end-neject+1:end) = [];
    t_cub_clean(:,i_nthreads) = col;
end

t_lin = t_lin_clean;
t_cub = t_cub_clean;

% Display the results
figure;
plot(number_of_threads_to_test, mean(t_lin), 'b-');
hold on;
plot(number_of_threads_to_test, mean(t_lin)-std(t_lin), 'b:');
plot(number_of_threads_to_test, mean(t_lin)+std(t_lin), 'b:');
xlabel('Number of threads');
ylabel('Average execution time');
title('Bilinear');

figure;
plot(number_of_threads_to_test, mean(t_cub), 'r*-');
hold on;
plot(number_of_threads_to_test, mean(t_cub)-std(t_cub), 'r:');
plot(number_of_threads_to_test, mean(t_cub)+std(t_cub), 'r:');
xlabel('Number of threads');
ylabel('Average execution time');
title('Bicubic');
