function [ P, Z, draw_L, NOISE, TTA_MU, TTA_SG, draw_PHI_MU, draw_PHI_SG, draw_PHI_COV, ...
    U_PHI_MU, U_PHI_SG, U_PHI_COV ] = model_variational( VARS, var_sizes, EPS1, EPS2 )

    % Unpack Variables
    [ NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, ...
        U_PHI_MU, U_PHI_SG, U_PHI_COV ] = unpack_vars( VARS, var_sizes );
    
    % Defaults
    if nargin < 3
        EPS1 = zeros(size(TTA_MU));
    end
    if nargin < 4
        EPS2_size = length(U_PHI_MU) + length(U_PHI_SG) + length(U_PHI_COV);
        EPS2 = zeros(1,EPS2_size);
    end
    
    % Sample from generator parameters
    idx_point = 1;
    % ---
    EPS2_1 = EPS2( idx_point:(idx_point+length(U_PHI_MU)-1) );
    draw_PHI_MU = PHI_MU + exp(U_PHI_MU).*EPS2_1;
    idx_point = idx_point + length(U_PHI_MU);
    % ---
    EPS2_2 = EPS2( idx_point:(idx_point+length(U_PHI_SG)-1) );
    draw_PHI_SG = PHI_SG + exp(U_PHI_SG).*EPS2_2;
    idx_point = idx_point + length(U_PHI_SG);
    % ---
    EPS2_3 = EPS2( idx_point:(idx_point+length(U_PHI_COV)-1) );
    draw_PHI_COV = PHI_COV + exp(U_PHI_COV).*EPS2_3;
    
    % Compute L
    N_params = var_sizes.s_TM(2);
    ind_offdiag = find(tril(ones(N_params,N_params),-1));
    draw_L = zeros(N_params,N_params);
    draw_L(ind_offdiag) = draw_PHI_COV;
    ind_diag = find(eye(N_params));
    draw_L(ind_diag) = exp(draw_PHI_SG);
    
    % Compute Z
    Z = TTA_MU + exp(TTA_SG).*EPS1;
    P = Z;
    
end

