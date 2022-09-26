function [ GRAD, grad_sizes, current_px, current_pth, current_pph, current_s_qph, current_g_qph, current_elbo ] = ...
    compute_gradients_eff( VARS, var_sizes, EPS1, EPS2, lambda, force_range, ALL_DATA)

    % RESHAPE GUIDE
    % P_vec = reshape(P, [ 1 N_params*N_subjects ]);
    % P = reshape(P_vec, [ N_subjects N_params ]);

    % Gradient Step
    ee = 0.0001;
    
    % Sizes
    N_subjects = var_sizes.s_TM(1);
    N_params = var_sizes.s_TM(2);
    
    % unpack NOISE
    [NOISE] = unpack_vars(VARS, var_sizes);
    N_noise = size(NOISE, 2);

    % P(TH|PHI) ----------------------------------------------------------
    
    [ ~, Z, L, ~, ~, ~, PHI_MU ] = ...
        model_variational( VARS, var_sizes, EPS1, EPS2 );
    base_pths = zeros(1, N_subjects);
    for j = 1:N_subjects
        base_pths(j) = log_pth(Z(j,:), L, PHI_MU, force_range);
    end
    base_pth = sum(base_pths);
    
    DPTH_DVARS = zeros(1, length(VARS));
    for i = 1:length(VARS)
        TEMP_VARS = VARS; TEMP_VARS(i) = TEMP_VARS(i) + ee;
        
        [ ~, this_Z, this_L, ~, ~, ~, this_PHI_MU] = ...
            model_variational( TEMP_VARS, var_sizes, EPS1, EPS2 );
        
        this_pths = zeros(1, N_subjects);
        for j = 1:N_subjects
            this_pths(j) = log_pth(this_Z(j,:), this_L, this_PHI_MU, force_range);
        end
        this_pth = sum(this_pths);
        this_g = (this_pth - base_pth)/ee;
        DPTH_DVARS(i) = this_g;
    end
    
    % P(PH) --------------------------------------------------------------
    
    [ ~, ~, L, ~, ~, ~, ~ ] = ...
        model_variational( VARS, var_sizes, EPS1, EPS2 );
    
    base_pph = log_pph(L, lambda);
    
    DPPH_DVARS = zeros(1, length(VARS));
    for i = 1:length(VARS)
        TEMP_VARS = VARS; TEMP_VARS(i) = TEMP_VARS(i) + ee;
        
        [ ~, ~, this_L, ~, ~, ~, ~] = ...
            model_variational( TEMP_VARS, var_sizes, EPS1, EPS2 );
        
        this_pph = log_pph(this_L, lambda);
        this_g = (this_pph - base_pph)/ee;
        DPPH_DVARS(i) = this_g;
    end
    
    % Q ------------------------------------------------------------------
    
    [~, ~, ~, ~, ~, TTA_SG, ~, ~, ~, U_PHI_MU, U_PHI_SG, U_PHI_COV ] = ...
        model_variational( VARS, var_sizes, EPS1, EPS2 );

    base_qphs = zeros(1, N_subjects);
    for j = 1:N_subjects
        base_qphs(j) = log_qph(EPS1(j,:),exp(TTA_SG(j,:))); % Subject
    end
    base_gen_qph = log_qph(EPS2, exp([U_PHI_MU, U_PHI_SG, U_PHI_COV])); % Generator
    base_qph = sum(base_qphs) + base_gen_qph;
    
    DQPH_DVARS = zeros(1, length(VARS));
    for i = 1:length(VARS)
        TEMP_VARS = VARS; TEMP_VARS(i) = TEMP_VARS(i) + ee;
        [~, ~, ~, ~, ~, this_TTA_SG, ~, ~, ~, this_U_PHI_MU, this_U_PHI_SG, this_U_PHI_COV ] = ...
            model_variational( TEMP_VARS, var_sizes, EPS1, EPS2 );
        this_qphs = zeros(1, N_subjects);
        for j = 1:N_subjects
            this_qphs(j) = log_qph(EPS1(j,:),exp(this_TTA_SG(j,:))); % Subject
        end
        this_gen_qph = log_qph(EPS2, exp([this_U_PHI_MU, this_U_PHI_SG, this_U_PHI_COV])); % Generator
        this_qph = sum(this_qphs) + this_gen_qph;

        this_g = (this_qph - base_qph)/ee;
        DQPH_DVARS(i) = this_g;
    end
    
    % PY -----------------------------------------------------------------
    
    P = model_variational( VARS, var_sizes, EPS1, EPS2 );
    
    % Get Current Py
    base_Px = zeros(N_subjects,1);
    parfor i = 1:N_subjects
        Data = ALL_DATA{i};
        base_Px(i) = log_py( P(i,:), Data, NOISE(i,:) );
    end
    
    % DPy/DNOISE
    DNOISE_vec = zeros(1, N_subjects*N_noise);
    parfor i = 0:(N_noise*N_subjects-1)
        
        row_idx = floor(i/N_noise) + 1;
        col_idx = mod(i,N_noise) + 1;
        
        TEMP_NOISE = NOISE(row_idx,:);
        TEMP_NOISE(col_idx) = TEMP_NOISE(col_idx) + ee;
        
        Data = ALL_DATA{row_idx};
        base_px = base_Px(row_idx);
        this_px = log_py(P(row_idx, :), Data, TEMP_NOISE );
        
        DNOISE_vec(i+1) = (this_px-base_px)/ee;
        
    end
    DNOISE = reshape(DNOISE_vec, [ N_noise N_subjects ])';
    DNOISE_DVARS = reshape(DNOISE, [ 1 N_noise*N_subjects ]);
    
    % Compute DPy/DP 
    DPX_vec = zeros(1,N_subjects*N_params);
    parfor i = 0:(N_params*N_subjects-1)
        
        row_idx = floor(i/N_params) + 1;
        col_idx = mod(i,N_params) + 1;
        
        TEMP_P = P;
        TEMP_P(row_idx, col_idx) = TEMP_P(row_idx, col_idx) + ee;
        
        Data = ALL_DATA{row_idx};
        base_px = base_Px(row_idx);
        this_px = log_py(TEMP_P(row_idx, :), Data, NOISE(row_idx,:));
        
        DPX_vec(i+1) = (this_px-base_px)/ee;
        
        % fprintf('R: %d, C: %d\n', row_idx, col_idx);
        
    end
   
    DPX = reshape(DPX_vec, [ N_params N_subjects ])';
    DPX_VEC = reshape(DPX, [ 1 N_params*N_subjects ]);
    
    % DP/DVARS
    base_P_vec = reshape(P, [ 1 N_params*N_subjects ]);
    
    DP_VEC = zeros(length(VARS), N_params*N_subjects);
    for i = 1:length(VARS)
        TEMP_VARS = VARS; TEMP_VARS(i) = TEMP_VARS(i) + ee;
        this_P = model_variational( TEMP_VARS, var_sizes, EPS1, EPS2 );
        this_P_vec = reshape(this_P, [ 1 N_params*N_subjects ]);
        this_g = (this_P_vec - base_P_vec)/ee;
        DP_VEC(i,:) = this_g;
    end
    
    % Compute DPy/DVARS
    DPX_DVARS = DPX_VEC*DP_VEC';
    
    % Store The noise gradient
    size_noise_variables = size(NOISE);
    n_noise_variables = size_noise_variables(1) * size_noise_variables(2);
    DPX_DVARS(1:n_noise_variables) = DNOISE_DVARS;
    
    % --------------------------------------------------------------------
    % Piece it all together:
    
    D_ELBO = DPX_DVARS + DPTH_DVARS + DPPH_DVARS - DQPH_DVARS;
    
    GRAD = D_ELBO';
    grad_sizes = var_sizes;
    
    current_px = sum(base_Px);
    current_pth = base_pth;
    current_s_qph = sum(base_qphs);
    current_g_qph = base_gen_qph;
    current_pph = base_pph;
    
    current_elbo = sum(base_Px) + base_pth + base_pph - base_qph;

end
