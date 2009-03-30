function [out] = Ave(z)
    out = (1/size(z,1))*sum(z);
end