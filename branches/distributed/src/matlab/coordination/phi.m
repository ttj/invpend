function [out] = phi(z, a, b)
    c = abs(a - b)/sqrt(4*a*b);

    out = (1/2)*(((a + b) * sigma_1(z + c)) + (a - b));
end