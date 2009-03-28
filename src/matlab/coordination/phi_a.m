function [out] = phi_a(z, r_sig, d_sig, h, a, b)
    out = rho_h(z / r_sig, h) * phi(z - d_sig, a, b);
end