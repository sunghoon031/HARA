clear all; close all; clc;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Create a synthetic dataset:
% (This is an example without the number of valid 2D-2D correspondences):

nViews = 100; 
connectivity = 0.2; 
outlier_ratio = 0.5;
inlier_noise_deg = 10;

[R_gt, edge_IDs, RR] = CreateSyntheticData(nViews, connectivity, outlier_ratio, inlier_noise_deg);



%% Run HARA:

% Set parameters:
median_thr = 1;
triplet_thr = 1;
err_sq_thr = 1;
nSamples = 10;
s_init = 10;
thrs_proportion = [0.1 0.2 0.3];
% nInliers_thrs = [5 0];

[R_est, time_initialization, time_optimization] = RunHARA(nViews, edge_IDs, RR, nSamples, s_init, thrs_proportion, triplet_thr, median_thr, err_sq_thr);
% [R_est, time_initialization, time_optimization] = RunHARA_usingNumberOfInlierMatches(nViews, edge_IDs, RR, nInliers_for_all_pairs, nSamples, s_init, thrs_proportion, triplet_thr, median_thr, err_sq_thr, nInliers_thrs)

disp(['Initialization took ', num2str(time_initialization), 's, Optimization took ', num2str(time_optimization), 's.'])  


%% Evaluate results:

[~,~, mean_error_L1_alignment, ~] = AlignRotationL1(R_gt, R_est);
[~, ~, rms_error_L2_alignment] = AlignRotationL2(R_gt, R_est);

disp(['Resulting errors (deg): theta_1 = ', num2str(mean_error_L1_alignment),...
    ', theta_2 = ', num2str(rms_error_L2_alignment)])





