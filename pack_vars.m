function [ VARS, sizes ] = pack_vars( NOISE, TTA_MU, TTA_SG, PHI_MU, PHI_SG, PHI_COV, ...
    U_PHI_MU, U_PHI_SG, U_PHI_COV )
    
    s_NO = size(NOISE);
    v_NO = reshape(NOISE,[s_NO(1)*s_NO(2) 1]);
    
    s_TM = size(TTA_MU);
    v_TM = reshape(TTA_MU,[s_TM(1)*s_TM(2) 1]);
    v_TS = reshape(TTA_SG,[s_TM(1)*s_TM(2) 1]);
    
    s_PM = size(PHI_MU);
    v_PM = reshape(PHI_MU,[s_PM(1)*s_PM(2) 1]);
    v_PS = reshape(PHI_SG,[s_PM(1)*s_PM(2) 1]);
    
    s_CV = size(PHI_COV);
    v_CV = reshape(PHI_COV, [s_CV(1)*s_CV(2) 1]);
    
    s_UPM = size(U_PHI_MU);
    v_UPM = reshape(U_PHI_MU, [s_UPM(1)*s_UPM(2) 1]);
    
    s_UPS = size(U_PHI_SG);
    v_UPS = reshape(U_PHI_SG, [s_UPS(1)*s_UPS(2) 1]);
    
    s_UCV = size(U_PHI_COV);
    v_UCV = reshape(U_PHI_COV, [s_UCV(1)*s_UCV(2) 1]);
    
    VARS = [v_NO; v_TM; v_TS; v_PM; v_PS; v_CV; v_UPM; v_UPS; v_UCV];
    
    sizes.s_TM = s_TM;
    sizes.s_PM = s_PM;
    sizes.s_NO = s_NO;
    sizes.s_CV = s_CV;
    sizes.s_UPM = s_UPM;
    sizes.s_UPS = s_UPS;
    sizes.s_UCV = s_UCV;
    
end

