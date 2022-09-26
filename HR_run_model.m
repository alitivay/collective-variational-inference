function [ Outputs ] = HR_run_model( Inputs, Parameters )

    Parameters = min(max(Parameters,0),1);
    
    % Time
    T_end = 180; % End time [minutes]
    Dt = 0.01; % Sampling time [minutes]
    time_range = 0:Dt:T_end;
    k_end = length(time_range);

    % Inputs
    UI = Inputs.Infusion.Values;
    UH = Inputs.Hemorrhage.Values;

    % Get initial conditions
    theta_s = HR_scale_parameters(Parameters);
    V0 = theta_s(11);
    H0 = theta_s(12);
    Va0 = 0.3*V0;
    Vv0 = 0.7*V0;
    Vr0 = H0*(Va0+Vv0);
    x0 = [Va0 Vv0 Vr0 0 0 0];

    % Simulation Init.
    x = x0;
    X = zeros(length(x0), k_end);
    Y = zeros(3, k_end);

    % Simulation Loop
    for k = 1:k_end
        t = time_range(k);
        X(:,k) = x;
        [ x, y ] = HR_transition_model( x, Parameters, [UI(k); UH(k)], Dt );
        Y(:,k) = y;
    end

    % Return Experiment Outputs
    Outputs.HCT = struct('Values', Y(1,:)', 'Times', time_range');
    Outputs.CO  = struct('Values', Y(2,:)', 'Times', time_range');
    Outputs.MAP = struct('Values', Y(3,:)', 'Times', time_range');
        
end

