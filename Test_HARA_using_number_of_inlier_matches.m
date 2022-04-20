clear all; close all; clc;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Load a real dataset:

dataset = 2;
% 1: Alamo, 2: Ellis_Island, 3: Gendarmenmarkt,  4: Madrid_Metropolis, 
% 5: Montreal_Notre_Dame, 6: Notre_Dame 1, 7: NYC_Library, 
% 8: Piazza_del_Popolo, 9: Piccadilly, 10: Roman_Forum, 
% 11: Tower_of_London, 12: Trafalgar, 13: Union_Square, 
% 14: Vienna_Cathedral, 15: Yorkminster,
    
[R_gt, edge_IDs, RR, nViews, nInliers_for_all_pairs] ...
    = LoadRealDataWithNumberOfInlierMatches(dataset);


%% Run HARA (using the number of inlier matches):

median_thr = 1;
triplet_thr = 1;
err_sq_thr = 1;
nSamples = 10;
s_init = 10;
thrs_proportion = [0.1 0.2 0.3];

nInliers_thrs = [5 0]; 


[R_est, time_initialization, time_optimization] = RunHARA_usingNumberOfInlierMatches(nViews, edge_IDs, RR, nInliers_for_all_pairs, nSamples, s_init, thrs_proportion, triplet_thr, median_thr, err_sq_thr, nInliers_thrs);

disp(['Initialization took ', num2str(time_initialization), 's, Optimization took ', num2str(time_optimization), 's.'])  


%% Evaluate results:

[~,~, mean_error_L1_alignment, ~] = AlignRotationL1(R_gt, R_est);
[~, ~, rms_error_L2_alignment] = AlignRotationL2(R_gt, R_est);

disp(['Resulting errors (deg): theta_1 = ', num2str(mean_error_L1_alignment),...
    ', theta_2 = ', num2str(rms_error_L2_alignment)])





