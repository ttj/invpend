function [out] = a_ij(q_i, q_j, r_sig, h, epsilon)
    out = rho_h( (sig_norm(q_i - q_j, epsilon)/r_sig), h);
    if out > 1 || out < 0
        'Error: adjacency matrix element out of range [0, 1]'
        out
    end
end