function [ lp ] = log_pth( z, L, mu_phi, force_range )
    
    % Generator-related density
    SIG = L*L';
    lpth = -(1/2)*(z-mu_phi)*inv(SIG)*(z-mu_phi)' - (1/2)*log(det(SIG));
    
    % Discourage z outside [0 1]
    slope = 100; scale = 100; lpc = 0;
    if force_range
        for zz = z
            if (zz <= (1/slope))
                lpc = lpc + slope*zz;
            elseif ( ((1/slope) < zz) && (zz < (1-(1/slope))) )
                lpc = lpc + 1;
            elseif (zz >= (1-(1/slope)))
                lpc = lpc + (-slope*zz + slope);
            end
        end
    end
    
    lp = lpth + scale * lpc;

end