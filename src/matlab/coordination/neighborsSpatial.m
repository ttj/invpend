function [Ni] = neighborsSpatial(i, q_i, q_js, r, d)
%Ni={j \in V : norm(q_j - q_i, 2) < r}
%  that is, the indices of j such that the distance between 
%  j and i is less than r
%
%  Recall that q_j, q_i \in R^m, so they are vectors containing the
%  positions

    kappa = r / d;

%Idea: take the overall state matrix A that contains all positions of all 
%      vertices (the q_j and q_i), then we have to do a search across every
%      element in this matrix to find those that are less than r apart
%      What search algorithm is best for this (noting that this will be 
%      computed many times)?

    idx = 0;
    Ni = zeros(size(q_js,1),1); %preallocate for largest possible size

    %iterate over all vertices (all rows in A)
    for j=1:size(q_js,1)
        x=norm(q_js(j,:) - q_i,2);
        if (x < r && j ~= i)
            %add to group
            idx = idx + 1;
            Ni(idx) = j;
        end
    end

    Ni = Ni(1:idx); %shrink idx
end