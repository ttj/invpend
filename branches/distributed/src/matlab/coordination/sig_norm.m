function [out] = sig_norm(z, epsilon)
    out = (1 / epsilon) * (sqrt(1 + epsilon * (norm(z,2))^2)-1);
    %out = norm(z,2);
end