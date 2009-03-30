function [out] = phi(z, a, b)
    %0 < a <= b
    c = abs(a - b)/sqrt(4*a*b);
    % => phi(0, a, b) = 0

    out = (1/2)*(((a + b) * sigma_1(z + c)) + (a - b));
end