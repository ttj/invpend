function [out] = sig_norm(z, epsilon)
    out = norm(z,2);
    %out = (1 / epsilon) * (sqrt(1 + epsilon * (norm(z,2))^2)-1);
end