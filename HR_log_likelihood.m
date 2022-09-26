function [ lpy ] = HR_log_likelihood( z, data, lognoise )
        
    lpy = 0;
    
    % Run the model
    Outputs = HR_run_model(data.Inputs, z);
    
    % Figure out the measured variables in this experiment
    measured_vars = fieldnames(data.Measurements);
    n_mv = numel(measured_vars);
    
    % Calculate log-likelihood w.r.t. available data in this experiment
    for k = 1:n_mv
        
        sigma = exp(lognoise(k));
        
        Yd = data.Measurements.(measured_vars{k}).Values;
        Td = data.Measurements.(measured_vars{k}).Times;
        
        Ym_raw = Outputs.(measured_vars{k}).Values;
        Tm_raw = Outputs.(measured_vars{k}).Times;
        Ym = interp1(Tm_raw,Ym_raw,Td);
        
        this_lpy = -(1/2)*sum((Yd-Ym).^2./sigma^2)-length(Yd)*log(sigma);
        lpy = lpy + this_lpy;
        
    end  

end

