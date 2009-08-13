function [q_val] = quantize(value, mode, precision)
    QUANT_MODE_FLOOR = 1;

    if mode == QUANT_MODE_FLOOR
        q_val = floor((2^(precision - 1)) * value)./(2^(precision - 1));
    else
        q_val = value;
    end
end