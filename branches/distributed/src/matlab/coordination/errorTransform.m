function [error] = errorTransform(framework, q, r_lattice, goal, leader)
    n = size(q,1); %# agents
    m = size(q,2); %# dimensions
    
    if framework == 2
        error = zeros(size(q));

        for i = 1 : n
            if i == leader
                %error(i) = 0;
                %error(i) = 0 + 0.01*sqrt(1 + norm(q(i,:) - goal, 2)) - 1;
                %error(i) = norm(q(i,:) - goal, 2);

                error(i) = q(i,:) - goal;
            elseif i < leader
                error(i) = q(i+1,:) - (q(i,:) + r_lattice);

                %error(i) = q(i,:) - (q(i-1,:) + r_lattice) + q(i,:) - goal -
                %i*r_lattice; (error converges to -r_lattice)
                %error(i) = q(i,:) - (q(i-1,:)) + q(i,:) - goal - i*r_lattice; %works: error converges to 0 and has a positive definite sum-of-squares Lyapunov function with a negative definite derivative
                %error(i) = 2*q(i,:) - q(i-1,:) - goal - i*r_lattice; %rewritten working version

                %can we rewrite goal and i*r_lattice in terms solely of the
                %other components?  I think so, let's see what might work to do
                %it

                %error to goal = current position - q(1,:) - goal + i*r_lattice
                %              = error(1) + i*r_lattice

                %error(i) = 2*q(i,:) - q(i-1,:) - goal - i*r_lattice; %this is basically the same as deviation energy

                %error(i) = q(i,:) - (q(i-1,:) + r_lattice) + 0.01*(sqrt(1 + norm(q(i,:) - goal,2))) - 1;

                %error(i) = q(i,:) - (q(i-1,:) + r_lattice) + norm(q(i,:) - goal, 2) - i*r_lattice;
            elseif i > leader
                error(i) = q(i,:) - (q(i-1,:) + r_lattice); %original transform
            end
        end

        %error(1) = error(1) + error(2);
        %error(1) = error(1) + sum(error(2:n));
    else
        error = zeros(n, m);
    end
end
