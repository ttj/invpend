function [ out ] = checkFlockingInvariant(q_last, q_now, spacing, error)
    N = size(q_now, 1);
    count = 0;
    count_group = N;
    all_last = 1;
    
    for i = 2 : N
        all_last = all_last && (norm(q_last(i) - q_last(i - 1) - spacing, 2) <= error);
    end
    
    for i = 2 : N
        %if last round was weak flock and this round is not, we have violated the invariant
        %also, this only satisfies for all nodes in a group, so it must be computed for each group, where a group is defined as
        %any N-hop connected group
        if all_last && (norm(q_last(i) - q_last(i - 1) - spacing, 2) <= error) && ~(norm(q_now(i) - q_now(i - 1) - spacing, 2) <= error)
            count = count + 1;
        end
    end

    out = (count == 0);
end
