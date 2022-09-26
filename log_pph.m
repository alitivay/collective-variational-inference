function [ logprob ] = log_pph( L, lambda )

        S = svd(L);
        logprob = - lambda * sum(S);

end

