function [epsilon] = proximityNetEdges(q, r, d)
  %epsilon(q) = {(i, j) \in V \times V : norm(q_j - q_i, 2) < r, i \neq j}
  
  epsilon = zeros(size(q,1)^2,size(q,2));

  size(q,1)
  for i=1:size(q,1)
      js = neighborsSpatial(i, q(i,:), q, r, d);
      for j=1:size(js)
        epsilon(i + j,:) = js';
      end
  end

end