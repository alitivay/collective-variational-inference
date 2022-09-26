
clear all;
close all;
clc;

% --- 1. Prepare and load data -------------------------------------------

DATASET = prepare_data;

% --- 2. Specify initial guesses based on the model ----------------------

% initial personalized model parameters
INIT.theta_mu_0 = [ 0.4018 0.1797 0.1802 0.1599 0.4560 0.3138 ...
    0.2603 0.3562 0.5616 0.1841 0.3921 0.3840 0.4488 0.2469 ];

% initial personalized model parameter uncertainty (logarithmic)
INIT.log_theta_sg_0 = -4;

% initial subject generator center
INIT.phi_mu_0 = INIT.theta_mu_0;

% initial subject generator spread (logarithmic)
INIT.log_phi_sg_0 = -1.6;

% initial signal quality parameters (noise power, logarithmic)
INIT.log_n_0 = [ -3.9 -1.6 2.3 ]; % [ sigma_HCT, sigma_CO, sigma_MAP ]

% --- 3. Specify inference options ---------------------------------------

% Maximum number of iterations
OPTIONS.n_iterations = 500;

% Adaptive moment estimation parameters (better not to change)
OPTIONS.alpha = 0.01;
OPTIONS.epsilon = 0.00000001;
OPTIONS.beta_1 = 0.5;
OPTIONS.beta_2 = 0.9;

% Subject generator formulation options
OPTIONS.diagonal_gaussian_generator = 0;
OPTIONS.uncertain_generator = 0;
OPTIONS.force_range = 1; % forces parameter range to [0 1]

% Compression parameter for the full-covariance gaussian subject generator
OPTIONS.lambda = 0;

% --- 4. Run collective variational inference ----------------------------

[RESULTS, HISTORY] = collective_vi(DATASET, INIT, OPTIONS);
save('RESULTS/CVI_RESULTS', 'RESULTS', 'HISTORY');
