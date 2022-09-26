function [ P, prob ] = model_generative( VARS, var_sizes, EPS )
    
    [NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV] = unpack_vars(VARS, var_sizes);
    
    % Compute L
    N_params = var_sizes.s_TM(2);
    ind_offdiag = find(tril(ones(N_params,N_params),-1));
    L = zeros(N_params,N_params);
    L(ind_offdiag) = PHI_COV;
    ind_diag = find(eye(N_params));
    L(ind_diag) = exp(PHI_SG);
    
    % Compute Z
    P = PHI_MU + EPS*L';
    
    % Compute Probability
    prob = exp(log_pth( P, L, PHI_MU, 1 ));
    
end

