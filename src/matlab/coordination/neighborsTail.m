function [Ni] = neighborsTail(i, q, r)
    if (i > size(q, 1) || i < 1) %index out of range
        Ni = [];
        return;
    end

    idx = 0;
    Ni = [];
    %Ni = zeros(size(q,1), 1); %preallocate for largest possible size
    
    %TODO: speed enhancement by [q_s, i_s] = sort(norm(q - q_i, 2)), check only nearest neighbors; might have to use map to do the norm like this
    
    %assume m = 1 (1-D)

    %iterate over all vertices
    for j = i + 1 : size(q,1)
        x = norm(q(j,:) - q(j - 1,:), 2);
        if ((x <= r) && (j ~= i))
            %add to group
            %idx = idx + 1;
            %Ni(idx) = j;
            Ni = j;
        else
            break; %stop on the first break in the chain of neighbors
        end
    end
    %Ni = Ni(1:idx); %shrink idx
end