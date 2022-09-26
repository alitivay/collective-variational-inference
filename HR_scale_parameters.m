function [ natural_params ] = hr_scale_parameters( scaled_params )

    % Parameter Bounds
    rlb = [ 0  0  0.01  10   1   0  5    0    0  0.001   1   0.1  3   70 ];
    rub = [ 6  6  0.5   100  10  2  100  500  3  0.05    5   0.5  6   130];
    
    natural_params = rlb + scaled_params .* (rub-rlb);
    
end

