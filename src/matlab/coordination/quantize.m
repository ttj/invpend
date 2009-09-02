function [u_new] = quantize(q_value, u_value, mode, precision)
    QUANT_MODE_FLOOR = 1;
    QUANT_MODE_BETA = 2;

    if mode == QUANT_MODE_FLOOR
        u_new = floor((2^(precision - 1)) * value)./(2^(precision - 1));
    elseif mode == QUANT_MODE_BETA
        if abs(u_value - q_value) < precision
            u_new = q_value; %don't change if smaller than quant_beta movement
        else
            u_new = u_value;
        end
    else
        u_new = value;
    end
end