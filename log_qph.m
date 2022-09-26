function [ lqph ] = log_qph( eps, sigma )
        
    lqph = -(1/2)*sum(eps.^2) - sum(log(sigma));

end

