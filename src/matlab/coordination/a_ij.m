function [out] = a_ij(q_i, q_j, r_sig, h, epsilon)
    out = rho_h( (sig_norm(q_i - q_j, epsilon)/r_sig), h);
end