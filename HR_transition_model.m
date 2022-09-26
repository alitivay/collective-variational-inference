function [ x_next, y ] = HR_transition_model( x, theta, u, Dt )
    
    % Read Inputs
    
    JI = u(1);
    JH = u(2);
    
    % Read States
    
    Va = x(:,1);
    Vv = x(:,2);
    Vr = x(:,3);
    rF = x(:,4);
    sR = x(:,5);
    sQ = x(:,6);
    
    % Read Parameters
    
    theta_s = HR_scale_parameters(theta);

    alpha_I = theta_s(:,1);
    alpha_H = theta_s(:,2);
    Kp      = theta_s(:,3);
    Kav     = theta_s(:,4);
    Kv      = theta_s(:,5);
    Kr      = theta_s(:,6);    
    tau_r   = theta_s(:,7);
    Kh      = theta_s(:,8);
    beta_v  = theta_s(:,9);
    KQ      = theta_s(:,10);
    V0      = theta_s(:,11);
    H0      = theta_s(:,12);
    Q0      = theta_s(:,13);
    Pa0     = theta_s(:,14);

    % Constants / Calculations
    
    Va0 = 0.3*V0;
    Vv0 = 0.7*V0;
    Pv0 = 6;
    
    H   = Vr./(Va+Vv);
    Ka  = Kv.*Kav;
    Rh_max = 5;
    
    % Pressures
    
    Pa = Pa0 + Ka .* (Va-Va0);
    Pv = Pv0 + Kv .* (Vv-Vv0);
    
    % Cardiac Output Model
    
    Q = Q0 + beta_v .* (Pv-Pv0) + sQ;
    sQ_next = sQ + Dt * ( -KQ.*(Q-Q0) );
    
    % Systemic Vascular Resistance Model
    
    R0 = (Pa0-Pv0)./Q0;
    R = R0 + Rh_max .* tanh((Kh./Rh_max).*(H-H0)) + sR;
    sR_next = sR + Dt .* ( -(1./tau_r).*sR - (Kr./tau_r).*(Pa-Pa0) );

    % Tissue Fluid Exchange Model
    
    JF = Kp .* (Va + Vv - Va0 - Vv0 - rF);
    rF_next = rF + Dt .* ( (1./(1+alpha_I)).*JI - (1./(1+alpha_H)).*JH );
    
    % Blood Circulation Model
    
    Va_next = Va + Dt .* ( Q - (Pa-Pv)./R - JH - JF );
    Vv_next = Vv + Dt .* ( - Q + (Pa-Pv)./R + JI );
    Vr_next = Vr + Dt .* ( - JH.*H );
    
    % Return Next States
    
    x_next = [ Va_next Vv_next Vr_next rF_next sR_next sQ_next ];
    
    y = [ H Q Pa ];
    
end
