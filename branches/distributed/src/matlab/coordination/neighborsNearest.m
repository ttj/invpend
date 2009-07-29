function [Ni] = neighborsNearest(i, q, r)
    if (i > size(q, 1)) %index out of range
        Ni = [];
        return;
    end

    idx = 1;
    Ni = zeros(size(q,1), 1); %preallocate for largest possible size
    
    %TODO: speed enhancement by [q_s, i_s] = sort(norm(q - q_i, 2)), check only nearest neighbors; might have to use map to do the norm like this
    
    %assume m = 1 (1-D)
    
    [q_s, i_s] = sort(q);

    %iterate over all vertices
    %for j = 1 : size(q,1)
    %    x = norm(q(j,:) - q(i,:), 2);
    %    if ((x <= r) && (j ~= i))
    %        %add to group
    %        idx = idx + 1;
    %        Ni(idx) = j;
    %    end
    %end
    
    if i >= 2 && (norm(q(i, :) - q(i - 1, :), 2) <= r)
        Ni(idx) = i - 1;
        idx = idx + 1;
    end
    
    if i <= size(q, 1) - 1 && (norm(q(i + 1, :) - q(i, :), 2) <= r)
        Ni(idx) = i + 1;
        idx = idx + 1;
    end

    Ni = Ni(1:idx-1); %shrink idx
end