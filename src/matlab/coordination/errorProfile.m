function [out] = errorProfile(q, e, d)
    N = size(q,1);
    
    e_factorL = ones(N,1) - 2.^(-([1:N]')); %1 - (1/(2^i))
    e_factorR = ones(N,1) - 2.^(-(N - ([1:N]'))); %1 - (1/(2^(N-i)))
    
    e = abs(e);
    d = abs(d);
    
    q_max = norm(q, inf);
    e_max = norm(e, inf);
    d_max = norm(d, inf);
    
    %idea: we need to normalize all of the errors based on the last nodes error from its final position, then we can apply this error profile
    %e_factorL
    %e_factorR
    %'e_profile'
    %e
    %d
    %e <= e_max.*e_factorL
    %e <= e_max.*e_factorR
    
    e;
    d;
    
    if (sum(sort(e) - e) ~= 0)
        'errors unsorted';

        if (sum(sort(d) - d) ~= 0)
            'not e => not d';
        end
    end

    if (sum(sort(d) - d) ~= 0)
        'deviation errors unsorted';
        
        if (sum(sort(e) - e) ~= 0)
            'not d => not e';
        end
    end
    
    %d <= d_max.*e_factor
end