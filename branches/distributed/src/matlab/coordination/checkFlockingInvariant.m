function [ out ] = checkFlockingInvariant(q_last, q_now, spacing, error)
    N = size(q_now, 1);
    count = 0;
    for i = 2 : N
        %if last round was weak flock and this round is not, we have violated the invariant
        if (norm(q_last(i) - q_last(i - 1) - spacing, 2) <= error) && ~(norm(q_now(i) - q_now(i - 1) - spacing, 2) <= error)
            count = count + 1;
        end
    end
    
    out = (count == 0);
end