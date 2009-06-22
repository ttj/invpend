function [error] = errorTransform(q, r_lattice, goal)
    n = size(q,1); %# agents
    m = size(q,2); %# dimensions
    error = zeros(size(q));

    for i = 1 : n
        if i == 1
            error(i) = q(i,:) - goal;
        else
            error(i) = q(i,:) - (q(i-1,:) + r_lattice);
        end
    end
end