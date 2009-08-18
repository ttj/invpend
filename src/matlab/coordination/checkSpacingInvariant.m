function [ out ] = checkSpacingInvariant(q, spacing)
    N = size(q, 1);
    count = 0;
    for i = 2 : N
        if (norm(q(i) - q(i - 1), 2) < spacing)
            count = count + 1;
            strcat('spacing: nodes ', int2str(i), ' and ', int2str(i-1), ' with spacing ', num2str(norm(q(i) - q(i - 1), 2)), '<', num2str(spacing))
        elseif q(i - 1) >= q(i)
            count = count + 1;
            strcat('collision: nodes ', int2str(i), ' and ', int2str(i-1), ' with spacing ', num2str(norm(q(i) - q(i - 1), 2)), '<', num2str(spacing))
        end
    end
    
    out = (count == 0);
end