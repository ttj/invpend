function E = deviationEnergy(framework, q, q_goal, r_comm, r_lattice, delta)
    n = size(q,1); %# agents
    m = size(q,2); %# dimensions
%     E = zeros(size(q));
%     for i=1:size(q,1)
%         %js=neighborsSpatialLattice(i, q(i,:), q, r, d, delta);
%         js=neighborsSpatial(i, q(i,:), q, r, d);
%         for j=1:size(js)
%             E(i,:) = E(i,:) + psi(norm(q(j,:) - q(i,:), 2) - d);
%         end
%     end
    %E = E .* (1/(norm(proximityNetEdges(q, r, d),1) + 1));

    %E = abs(q(:,:) - q_goal(:,:));
    fv = (0:1:n-1)'*r_lattice
    %E = abs(q(:,:) - q(1,:) - fv)
    E = q(:,:) - q(1,:) - fv
    E(1) = 0;

    %E = q(:,:) - q_goal(:,:);
    
    if framework == 3
        E = q(:,:);
    end
end
