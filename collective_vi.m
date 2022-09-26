function [ RESULTS, HISTORY ] = collective_vi(DATASET, INIT, OPTIONS)

    % Parameter Initializations
    theta_mu_0 = INIT.theta_mu_0;
    phi_mu_0 = INIT.phi_mu_0;
    log_theta_sg_0 = INIT.log_theta_sg_0;
    log_phi_sg_0 = INIT.log_phi_sg_0;
    log_n_0 = INIT.log_n_0;

    % Optimization Parameters
    n_iterations = OPTIONS.n_iterations;
    lambda = OPTIONS.lambda;
    alpha = OPTIONS.alpha;
    epsilon = OPTIONS.epsilon;
    beta_1 = OPTIONS.beta_1;
    beta_2 = OPTIONS.beta_2;
    diagonal_gaussian_generator = OPTIONS.diagonal_gaussian_generator;
    uncertain_generator = OPTIONS.uncertain_generator;
    force_range = OPTIONS.force_range;

    % Start parallel cluster for higher speed
    if isempty(gcp('nocreate'))
        this_cluster = parcluster('local');
        this_cluster.NumWorkers = 6;
        this_pool = parpool(this_cluster, this_cluster.NumWorkers);
    end

    % Calculate Vector Sizes
    N_params = length(theta_mu_0);
    N_subjects = length(DATASET);
    N_cov = (N_params^2-N_params)/2; % Generator Covariance Parameter Size
    N_gen_params = N_params*2+N_cov; % All Generator Parameter Size

    % Initialize Parameters
    TTA_MU = repmat(theta_mu_0,N_subjects,1); % Subject-Specific Posterior Means
    TTA_SG = repmat(log_theta_sg_0*ones(1,N_params),N_subjects,1); % Subject-Specific Posterior STDs
    PHI_MU = repmat(phi_mu_0,1,1); % Population Posterior Means
    PHI_SG = repmat(log_phi_sg_0*ones(1,N_params),1,1); % Population Posterior STDs
    NOISE  = repmat(log_n_0, N_subjects, 1);
    PHI_COV = zeros(1,N_cov); % Generator Covariance: Parameter Vector
    U_PHI_MU  = -3*ones(1,N_params);
    U_PHI_SG  = -3*ones(1,N_params);
    U_PHI_COV = -3*ones(1,N_cov);

    % History Parameters
    H_TTA_MU = zeros([size(TTA_MU) n_iterations]);   % Subject-Specific Posterior Means
    H_TTA_SG = zeros([size(TTA_SG) n_iterations]);   % Subject-Specific Posterior STDs
    H_PHI_MU = zeros([size(PHI_MU) n_iterations]);   % Population Means
    H_PHI_SG = zeros([size(PHI_SG) n_iterations]);   % Population STDs
    H_PHI_COV = zeros([size(PHI_COV) n_iterations]); % Population Covariances
    H_ELBO = zeros([N_subjects n_iterations]); % Subject-Separated ELBO Values
    H_NOISE  = zeros([size(NOISE), n_iterations]); % Logarithm of Noise Sigma

    % Pack Optimization Variables
    [VARS, ~] = pack_vars(NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, U_PHI_MU, U_PHI_SG, U_PHI_COV);

    % Initialize Moments
    M = zeros(length(VARS),1);
    V = zeros(length(VARS),1);
    T = 0;

    % Initialize Progress Plot
    figure;
    set(gcf, 'Position', [100, 100, 600, 350]);
    xlabel('Iteration'); ylabel('ELBO'); grid on;

    for idx = 1:n_iterations

        T = T+1;

        % Stochasticity: Pick Random Epsilon (for model parameters)
        EPS1 = randn(length(DATASET), N_params);
        
        % Stochasticity: Pick Random Epsilon (for uncertain generator)
        if uncertain_generator
            EPS2 = randn(1, N_gen_params);
        else
            EPS2 = zeros(1, N_gen_params);
        end

        % Get Current Gradient
        ta = tic;
        [VARS, var_sizes] = ...
            pack_vars(NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, U_PHI_MU, U_PHI_SG, U_PHI_COV);
        [GRAD, ~, current_px, current_pth, current_pph, current_s_qph, current_g_qph, current_elbo ] = ...
            compute_gradients_eff(VARS, var_sizes, EPS1, EPS2, lambda, force_range, DATASET);
        tb = toc(ta);

        % Display Stats
        fprintf('Iteration: %d, ELBO: %f, PY: %f, PTH: %f, PPH: %f, S-QPH: %f, G-QPH: %f\n', ...
            idx, current_elbo, current_px, current_pth, current_pph, current_s_qph, current_g_qph);
        fprintf('Time: %f\n', tb);

        % Moments
        M = beta_1*M + (1-beta_1)*GRAD;
        V = beta_2*V + (1-beta_2)*GRAD.^2;

        % Bias-Correction
        M_hat = M/(1-beta_1^T);
        V_hat = V/(1-beta_2^T);

        % Gradient Step
        UPDATE = M_hat./(sqrt(V_hat)+epsilon);
        VARS = VARS + alpha * UPDATE;

        % Unpack Optimization Variables
        [ NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, U_PHI_MU, U_PHI_SG, U_PHI_COV ] = unpack_vars( VARS, var_sizes );

        % Force Diagonal Gaussian Generator
        if diagonal_gaussian_generator
            PHI_COV = zeros(1,N_cov);
            U_PHI_COV = -10*ones(1,N_cov);
        end
        
        % Force Certain Generator
        if ~uncertain_generator
             U_PHI_MU  = -10*ones(1,N_params);
             U_PHI_SG  = -10*ones(1,N_params);
             U_PHI_COV = -10*ones(1,N_cov);
        end

        % Save Optimization Path
        H_TTA_MU(:,:,idx) = TTA_MU;
        H_TTA_SG(:,:,idx) = TTA_SG;
        H_PHI_MU(:,:,idx) = PHI_MU;
        H_PHI_SG(:,:,idx) = PHI_SG;
        H_PHI_COV(:,:,idx) = PHI_COV;
        H_NOISE(:,:,idx) = NOISE;
        H_ELBO(:,idx) = current_elbo;

        if mod(idx,1)==0
            plot(idx,mean(current_elbo),'ob'); hold on;
            xlabel('Iteration'); ylabel('ELBO'); grid on;
            drawnow();
        end

    end
    
    % Return Outputs
    RESULTS.TTA_MU  = TTA_MU;
    RESULTS.TTA_SG  = TTA_SG;
    RESULTS.PHI_MU  = PHI_MU;
    RESULTS.PHI_SG  = PHI_SG;
    RESULTS.PHI_COV = PHI_COV;
    RESULTS.NOISE   = NOISE;
    RESULTS.VARS    = VARS;
    RESULTS.var_sizes = var_sizes;
    
    % Return History
    HISTORY.H_TTA_MU  = H_TTA_MU;
    HISTORY.H_TTA_SG  = H_TTA_SG;
    HISTORY.H_PHI_MU  = H_PHI_MU;
    HISTORY.H_PHI_SG  = H_PHI_SG;
    HISTORY.H_PHI_COV = H_PHI_COV;
    HISTORY.H_NOISE   = H_NOISE;
    HISTORY.ELBO      = H_ELBO;

end
