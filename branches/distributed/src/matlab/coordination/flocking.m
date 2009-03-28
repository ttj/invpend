function [ out ] = flocking(N, m, coord_min, coord_max, r, d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %original sim:
    %d=7
    %r=1.2*d=8.4
    %d'=0.6*d=
    %r'=1.2*d'=
    %epsilon=0.1
    %a=b=5 for \phi
    %h=0.2 for \phi_a
    %h=0.9 for \phi_b
    %step size=0.01 to 0.03s (33hz-100hz)
    %N=150
    %initial positions \in gaussian with var=2500
    %initial velocities \in [-2, -1]^2 box
    
    epsilon = 0.1;
    a=5;
    b=5;
    ha=0.2;
    hb=0.9;
    
    Ts=0.01;
    Tc=0.01;
    fs=1/Ts;
    fc=1/Tc;

    kappa = r / d;
    
    vel_min = -2;
    vel_max = -1;
    %generate velocity matrix
    p = vel_max + (vel_min - vel_max).*rand(N, m);

    %generate state matrix randomly with N nodes \in R^m
    %q = coord_max + (coord_min - coord_max).*rand(N, m);
    
    %generate state matrix randomly such that no nodes are already
    %neighbors
    q=zeros(N, m);
    for i=1:1:N
        q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);
        
        for j=1:1:i %only need to go up to i's already initialized
            if i ~= j
                %ensure no vertices start as neighbors, that is 
                %(norm(q_i,q_j,2) < r) != true
                while size(neighborsSpatial(i, q(i,:), q, r, d),1) > 0
                    q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);
                end
            end
        end
    end

    spatial_neighbors = zeros(N, N, N);

    tcyc=Tc;
    tdiv=5; %minimum is 1 (1 division per control cycle, same as only looking at control cycle)
    tmax=0.5;
    time_ctrl=[0:tcyc:tmax]';      %vector of control cycle times
    time_traj=[0:tcyc/tdiv:tmax]'; %vector of trajectories
    u_steps=(tmax-tcyc)/tcyc+1;     %number of control update steps
    steps=tmax/(tcyc/tdiv); %total numbers of trajectory steps
    updates=10; %number of times to update the plots
    
    %system evolution
    for t=0:tcyc:tmax-tcyc
        t_i=round(t/tcyc)+1;
        
        if mod(t_i-1, u_steps / updates) == 0
            %we will call the below plotting loop every some odd iterations of the
            %system evolution to see how it's doing
            if m == 2 || m == 3
                figure;
                hold on;
            end

            if m == 2
                scatter(q(:,1),q(:,2));
            elseif m == 3
                scatter3(q(:,1),q(:,2),q(:,3));
            end

            %generate neighbors for all nodes and plot
            %TODO: have this simply create a data structure and then plot that
            %      What will that data structure look like?  Rather complicated
            %      multi-dimensional object without consistent number of elements
            %TODO: move to function
            for i=1:size(q,1)
                neighbors_i = neighborsSpatial(i, q(i,:), q, r, d); %exclude self here if desired q_js: cat(1,q(1:i,:),q(i+1:size(q,1),:))
                %spatial_neighbors(:,:,i) = neighbors_i;

                %TODO: speed this up by removing the for loop and doing
                %everything via index operations after generating the neighbors
                for j=1:size(neighbors_i)
                    if m == 2
                        plot([q(i,1) q(neighbors_i(j,1),1)], [q(i,2) q(neighbors_i(j,1),2)]);
                    elseif m == 3
                        plot3([q(i,1) q(neighbors_i(j,1),1)], [q(i,2) q(neighbors_i(j,1),2)], [q(i,3) q(neighbors_i(j,1),3)]);
                    end
                end
            end

            %spatial_neighbors
        end
        
        %compute control (based on state vector, possibly delayed, etc)
        %TODO: add proper control, random probably isn't effective ;-)
        u = (vel_max + (vel_min - vel_max).*rand(N, m))./10;
        %u = zeros(N, 1);
        
        %u_i = u_i^\alpha + u_i^\gamma
        %n_ij = ((q_j - q_i ) / (sqrt(1 + epsilon * norm(q_j - q_i, 2))^2));
        
        %rho_h(z) = 1                                    if z \in [0,h)
        %           1/2 * (1 + cos(pi * (z - h)/1 -h))   if z \in [h,1]
        %           0                                    else
        %for h \in (0,1)
        
        %a_ij(q) = rho_h( sig_norm ( q_j - q_i ) / r_sig )
        %and a_ii(q) = 0 for all i, q
        %a_ij(q) \in [0, 1]
        % take h =1 (doesn't this violate set definition?
        % r_sig = sig_norm( r )
        r_sig = sig_norm(r, epsilon);
        d_sig = sig_norm(d, epsilon);
        
        %sig_norm(z) = (1 / epsilon) * (sqrt(1 + epsilon * (norm(z,2))^2)-1)
        
        %u_i = sumNghbs(\phi_a(sig_norm(q_j-q_i)) * n_ij) + 
        %      sumNghbs(a_ij(q) * (p_j - p_i))
        
        for i=1:N
            js = neighborsSpatial(i, q(i,:), q, r, d);
            for j=1:size(js,1)
                %TODO: verify equations and all functions
                u(i) = phi_a(sig_norm ( q(j,:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * ((norm(q(j,:) - q(i,:),2) ) / (sqrt(1 + epsilon * norm(q(j,:) - q(i,:), 2))^2));
                u(i) = u(i) + rho_h( sig_norm ( q(j,:) - q(i,:), epsilon ) / r_sig, ha) * norm((p(j,:) - p(i,:)),2);
                %TODO: add goal term (u_gamma)
            end
        end

        for t_j=(t_i-1)*tdiv+2 : 1 : tdiv+(t_i-1)*tdiv+1
            tt=t_j*(tcyc/tdiv);
            
            %run system evolution
            %use matrix exp?
            
            %q'=p;     => [q'; p']=[0 1; 0 0]*[q; p] + [0;1]*[u]
            %p'=u;
            %TODO: this is currently wrong, but it does result in some
            %notion of stability (convolution integral?)
            q=exp(p*tt).*q;
            p=exp(p*tt).*u;
        end
    end
end
