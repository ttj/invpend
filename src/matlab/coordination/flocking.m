function [ out ] = flocking(N, m, coord_min, coord_max, r, d, tdiv, tmax, updates, plotControls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%       tdiv:       minimum is 1 (1 division per control cycle, same as 
%                   only looking at control cycle)
%       tmax:       time to run the system for (t_final)
%       updates:    number of times to update the plots
%       plotControls:   if 1, plots all controls as functions of time

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
    
    c1=1.75;
    c2=0.5;
    
    epsilon = 0.1;
    a=5;
    b=5;
    ha=0.2;
    hb=0.9;
    delta=d/10;

    r_sig = sig_norm(r, epsilon);
    d_sig = sig_norm(d, epsilon);
    
    Ts=0.01;
    Tc=0.01;
    fs=1/Ts;
    fc=1/Tc;

    kappa = r / d;
    
    vel_min = -2;
    vel_max = -1;

    %generate velocity matrix
    p = vel_max + (vel_min - vel_max).*rand(N, m);
    
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
    time_ctrl=[0:tcyc:tmax-tcyc]';       %vector of control cycle times
    time_traj=[0:tcyc/tdiv:tmax]';  %vector of trajectories
    u_steps=((tmax-tcyc)/tcyc)+1;   %number of control update steps
    steps=tmax/(tcyc/tdiv);         %total numbers of trajectory steps

    %goals for gamma control term
    qd=zeros(N, m);         %preallocate
    pd=ones(N, m).*(d/10);
    
    %initial control value
    %u = (vel_max + (vel_min - vel_max).*rand(N, m));
    u = zeros(N, m);

    %generate an equally spaced (by d) rectangular grid goal
    %TODO: generalize for m-dimensional (use cat and m-th root of N [need m loops?])
    for i=1:ceil(sqrt(N))
        for j=1:ceil(sqrt(N))
            qd((i-1)*ceil(sqrt(N))+j,:) = [((i-1)*d) ((j-1)*d)] + [coord_max coord_max];
        end
    end
    qd=qd(1:N,:); %shrink to maximum N
    scatter(qd(:,1),qd(:,2));
    
    qr=qd;
    pr=pd;
    
    %start at goal
    %q=qd;
    %p=pd;
    
    subplotRows=1;%ceil(sqrt(updates));
    subplotCols=1;%floor(sqrt(updates));
    subplotCount=1;

    u_history = zeros([round(u_steps), N, m]);
    q_history = zeros([round(steps), N, m]);
    p_history = zeros([round(steps), N, m]);
    qr_history = zeros([round(steps), N, m]);
    pr_history = zeros([round(steps), N, m]);
    
    %system evolution
    for t=0:tcyc:tmax-tcyc
        t_i=round(t/tcyc)+1;
        if t_i==0 || mod(t_i+1, u_steps / updates) == 0
            %we will call the below plotting loop every some odd iterations of the
            %system evolution to see how it's doing
            if m == 2 || m == 3
                figure;
                hold on;
            end

            if m == 2
                subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),q(:,2));
                subplot(subplotRows,subplotCols,subplotCount), scatter(qd(:,1)+pd(:,1).*t,qd(:,2)+pd(:,2).*t,'r');
                %out=[[qd(:,1) (qd(:,1)+pd(:,1))] [qd(:,2) (qd(:,2)+pd(:,2))]];
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplotCount = subplotCount + 1;
            elseif m == 3
                subplot(subplotRows,subplotCols,subplotCount), scatter3(q(:,1),q(:,2),q(:,3));
                subplot(subplotRows,subplotCols,subplotCount), scatter3(qd(:,1),qd(:,2),qd(:,3),'r');
                %subplotCount = subplotCount + 1;
            end

            %generate neighbors for all nodes and plot
            %TODO: have this simply create a data structure and then plot that
            %      What will that data structure look like?  Rather complicated
            %      multi-dimensional object without consistent number of elements
            %TODO: move to function
            for i=1:size(q,1)
                neighbors_i = neighborsSpatialLattice(i, q(i,:), q, r, d, delta); %exclude self here if desired q_js: cat(1,q(1:i,:),q(i+1:size(q,1),:))
                %spatial_neighbors(:,:,i) = neighbors_i;

                %TODO: speed this up by removing the for loop and doing
                %everything via index operations after generating the neighbors
                for j=1:size(neighbors_i)
                    %prune data
                    %if j > i
                    %    continue;
                    %end
                    
                    if m == 2
                        plot([q(i,1) q(neighbors_i(j,1),1)], [q(i,2) q(neighbors_i(j,1),2)]);
                    elseif m == 3
                        plot3([q(i,1) q(neighbors_i(j,1),1)], [q(i,2) q(neighbors_i(j,1),2)], [q(i,3) q(neighbors_i(j,1),3)]);
                    end
                end
            end

            %spatial_neighbors
        end
        
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
        
        %sig_norm(z) = (1 / epsilon) * (sqrt(1 + epsilon * (norm(z,2))^2)-1)
        
        %u_i = sumNghbs(\phi_a(sig_norm(q_j-q_i)) * n_ij) + 
        %      sumNghbs(a_ij(q) * (p_j - p_i))
        
        %store all controls over time
        u_history(t_i,:,:) = u(:,:);
        
        %reinitialze all controls (not dependent upon past control value)
        u = zeros(N, m);
        
        %compute control (based on state vector, possibly delayed, etc)
        for i=1:N
            js = neighborsSpatial(i, q(i,:), q, r, d);
%            if size(js,1) > -1
                for j=1:size(js,1)
                    %TODO: verify equations and all functions
                    u(i,:) = u(i,:) + phi_a(sig_norm ( q(j,:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * (((q(j,:) - q(i,:)) ) / (sqrt(1 + epsilon * ((norm(q(j,:) - q(i,:), 2))^2))));
                    u(i,:) = u(i,:) + (a_ij(q(i,:), q(j,:), r_sig, ha, epsilon) * ((p(j,:) - p(i,:))));
                end
%             else
%                 %no neighbors yet: have to compute something
%                 %TODO: fix, not right (always 0 obviously)
%                 u(i,:) = u(i,:) + phi_a(sig_norm ( q(i,:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * (((q(i,:) - q(i,:)) ) / (sqrt(1 + epsilon * norm(q(i,:) - q(i,:), 2))^2));
%                 u(i,:) = u(i,:) + (a_ij(q(i,:), q(i,:), r_sig, ha, epsilon) * ((p(i,:) - p(i,:))));
%             end
            
            %add gamma goal term
            u(i,:) = u(i,:) - c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));
        end

        for t_j=(t_i-1)*tdiv+1 : 1 : tdiv+(t_i-1)*tdiv+1
            tt=t_j*(tcyc/tdiv);
            
            %run system evolution
            
            %state space form
            %q'=p;     => [q'; p']=[0 1; 0 0]*[q; p] + [0;1]*[u]
            %p'=u;
            
            %solution is:
            %  expm(A*(t-t0))*x0+int(expm(A*(t-t0))*B*u,tau,t0,t)
            %for
            %  syms q p q0 p0 tau t t0
            %  x = [q;p]
            %  x0 = [q0; p0]
            %  A = [0 1; 0 0]
            %  B = [0; 1]
            %which expands to
            %  [x]=[q;= [q0 + p0*(t - t0) + u*(t - t0)^2;
            %       p]   p0 + u*(t - t0)]
            
            %store all state variables over time
            q_history(t_j,:,:) = q(:,:);
            p_history(t_j,:,:) = p(:,:);
            qr_history(t_j,:,:) = qr(:,:);
            pr_history(t_j,:,:) = pr(:,:);
            
            %TODO: verify correctness of this as a solution
            %      in general, can we use this solution if the control
            %      is nonlinear?
            q=q + p.*(tcyc/tdiv) + u.*((tcyc/tdiv)^2);
            p=p + u.*(tcyc/tdiv);
            
            %gamma agent
            qr=qr + pr.*(tcyc/tdiv) + fr(qr, pr).*((tcyc/tdiv)^2);
            pr=pr + fr(qr, pr).*(tcyc/tdiv);
        end
    end
    
    out = u_history;

    %plot controls over time
    if plotControls == 1
        for i=1:N
            figure;
            hold on;
            plot(time_ctrl,u_history(:,i,1),'b');
            plot(time_ctrl,u_history(:,i,2),'r');
            plot(time_ctrl,sqrt(u_history(:,i,1).^2 + u_history(:,i,2).^2),'g');
            figure;
            hold on;
            plot(time_traj,q_history(:,i,1) - qd(i,1),'b');
            plot(time_traj,q_history(:,i,2) - qd(i,1),'r');
            plot(time_traj,p_history(:,i,1) - pd(i,1),'g');
            plot(time_traj,p_history(:,i,2) - pd(i,1),'c');
        end
    end

    %average position and velocity for moving reference frame
    qc=Ave(q)
    pc=Ave(p)
end


%1: fix lattice: goals independently in grid with multiple objectives

%2: 2 particle set up with one goal: how are these aligned? simplify
%functions enough for proofs

%3: obstacles

%4: failure with stuck value, starting in lattice

%5: one node slowly updating, everyone else ignores eventually and treats
%it as an obstacle: growing larger and larger obstacle as uncertainty is
%compounded over and over again

%john sicyclis: most general necessary and sufficient: on asynchronous
%iterated systems

%mani chandy's paper was like this going towards message passing systems
