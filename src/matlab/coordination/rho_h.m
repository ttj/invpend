function [out] = rho_h(z, h)
    %rho_h(z) = 1                                    if z \in [0,h)
    %           1/2 * (1 + cos(pi * (z - h)/1 -h))   if z \in [h,1]
    %           0                                    else
    %for h \in (0,1)
    if (z >= 0 && z < h)
        out = 1;
    elseif (z >= h && z <= 1)
        out = (1/2) * (1 + cos(pi * (z - h)/1 -h));
    else
        out = 0;
    end
end