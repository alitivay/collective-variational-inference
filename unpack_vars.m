function [ NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, ...
    U_PHI_MU, U_PHI_SG, U_PHI_COV ] = unpack_vars( VARS, sizes )
    
    idx_point = 1;
    
    % Unpack Noise
    length_NO = sizes.s_NO(1)*sizes.s_NO(2);
    v_NO = VARS( idx_point:(idx_point+length_NO-1) );
    idx_point = idx_point + length_NO;
    
    % Unpack Theta Parameters
    length_TM = sizes.s_TM(1)*sizes.s_TM(2);
    v_TM = VARS( idx_point:(idx_point+length_TM-1) );
    idx_point = idx_point + length_TM;
    % ---
    v_TS = VARS( idx_point:(idx_point+length_TM-1) );
    idx_point = idx_point + length_TM;
    
    % Unpack Phi Parameters
    length_PM = sizes.s_PM(1)*sizes.s_PM(2);
    v_PM = VARS( idx_point:(idx_point+length_PM-1) );
    idx_point = idx_point + length_PM;
    % ---
    v_PS = VARS( idx_point:(idx_point+length_PM-1) );
    idx_point = idx_point + length_PM;
   
    % Unpack Covariance Parameters
    length_CV = sizes.s_CV(1)*sizes.s_CV(2);
    v_CV = VARS( idx_point:(idx_point+length_CV-1) );  
    idx_point = idx_point + length_CV;
    
    % Unpack Generator Uncertainties:
    length_UPM = sizes.s_UPM(1)*sizes.s_UPM(2);
    v_UPM = VARS( idx_point:(idx_point+length_UPM-1) ); 
    idx_point = idx_point + length_UPM;
    % ---
    length_UPS = sizes.s_UPS(1)*sizes.s_UPS(2);
    v_UPS = VARS( idx_point:(idx_point+length_UPS-1) ); 
    idx_point = idx_point + length_UPS;
    % ---
    length_UCV = sizes.s_UCV(1)*sizes.s_UCV(2);
    v_UCV = VARS( idx_point:(idx_point+length_UCV-1) ); 
    
    % Reshaped Outputs
    NOISE = reshape(v_NO, sizes.s_NO);
    TTA_MU = reshape(v_TM, sizes.s_TM);
    TTA_SG = reshape(v_TS, sizes.s_TM);
    PHI_MU = reshape(v_PM, sizes.s_PM);
    PHI_SG = reshape(v_PS, sizes.s_PM);
    PHI_COV = reshape(v_CV, sizes.s_CV);
    U_PHI_MU = reshape(v_UPM, sizes.s_UPM);
    U_PHI_SG = reshape(v_UPS, sizes.s_UPS);
    U_PHI_COV = reshape(v_UCV, sizes.s_UCV);
    
end

