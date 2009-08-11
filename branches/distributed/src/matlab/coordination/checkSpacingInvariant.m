function [ out ] = checkSpacingInvariant(q, spacing)
    N = size(q, 1);
    count = 0;
    for i = 2 : N
        if (norm(q(i) - q(i - 1), 2) < spacing) || q(i - 1) >= q(i)
            count = count + 1;
        end
    end
    
    out = (count == 0);
end