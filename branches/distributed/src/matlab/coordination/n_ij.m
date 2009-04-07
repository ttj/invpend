function [out] = n_ij(q_i, q_j, epsilon)
    out = (q_j - q_i) / sqrt(1 + (epsilon * ((norm(q_j - q_i, 2))^2)));
end