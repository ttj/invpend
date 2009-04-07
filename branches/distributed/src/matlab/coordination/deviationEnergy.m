function E = deviationEnergy(q, r, d, delta)
    %n = size(q,1); %# agents
    %m = size(q,2); %# dimensions
    E = zeros(size(q));
    for i=1:size(q,1)
        %js=neighborsSpatialLattice(i, q(i,:), q, r, d, delta);
        js=neighborsSpatial(i, q(i,:), q, r, d);
        for j=1:size(js)
            E(i,:) = E(i,:) + psi(norm(q(j,:) - q(i,:), 2) - d);
        end
    end
    %E = E .* (1/(norm(proximityNetEdges(q, r, d),1) + 1));
end