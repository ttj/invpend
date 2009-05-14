function [out] = sig_norm_old(z, epsilon)
    out = (1 / epsilon) * (sqrt(1 + epsilon * (norm(z,2))^2)-1);
end