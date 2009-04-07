function [ out ] = flocking(N, m, coord_min, coord_max, r, d, tdiv, tmax, updates, plotControls)
    %close all;
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
    
    %c1=1.75;
    %c2=25;
    c1=0.25;
    c2=0.1;
    
    epsilon = 0.1;
    a=5;
    b=5;
    ha=0.2;
    hb=0.9;
    delta=d/5;

    r_sig = sig_norm(r, epsilon);
    d_sig = sig_norm(d, epsilon);
    
    Ts=0.01;
    Tc=0.01;
    fs=1/Ts;
    fc=1/Tc;

    kappa = r / d;

    vel_min = 10;
    vel_max = 0;

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
    pd=ones(N, m);
    
    %initial control value
    %u = (vel_max + (vel_min - vel_max).*rand(N, m));
    u = zeros(N, m);
    uGradient = zeros(N, m);
    uConsensus = zeros(N, m);
    uGamma = zeros(N, m);

    %generate an equally spaced (by d) rectangular grid goal
    %TODO: generalize for m-dimensional (use cat and m-th root of N [need m loops?])
    if m == 1
        for i=1:N
            qd(i,:) = ((i-1)*d) + coord_max;
        end
        scatter(qd(:,1),zeros(N,1),[1:N]');
    elseif m == 2
        for i=1:ceil(sqrt(N))
            for j=1:ceil(sqrt(N))
                qd((i-1)*ceil(sqrt(N))+j,:) = [((i-1)*d) ((j-1)*d)] + [coord_max coord_max];
            end
        end
        scatter(qd(:,1),qd(:,2));
    elseif m == 3
        for i=1:ceil(N^(1/3))
            for j=1:ceil(N^(1/3))
                for k=1:ceil(N^(1/3))
                    qd((i-1)*ceil(N^(2/3))+(j-1)*ceil(N^(1/3))+k,:) = [((i-1)*d) ((j-1)*d) ((k-1)*d)] + [coord_max coord_max coord_max];
                end
            end
        end
        scatter(qd(:,1),qd(:,2),qd(:,3));
    end
    qd=qd(1:N,:); %shrink to maximum N
    %qd=qd.*5;
    qd=ones(N,m).*coord_max*5;
    %pd=zeros(N,m);
    
    qr=qd;
    pr=pd;
    
    %start at goal
    %q=qd;
    %p=pd;
    
    de = deviationEnergy(q,r,d,delta)
    %pne = proximityNetEdges(q,r,d)
    subplotRows=1;%ceil(sqrt(updates));
    subplotCols=1;%floor(sqrt(updates));
    subplotCount=1;

    u_history = zeros([round(u_steps), N, m]);
    uGradient_history = zeros([round(u_steps), N, m]);
    uConsensus_history = zeros([round(u_steps), N, m]);
    uGamma_history = zeros([round(u_steps), N, m]);
    q_history = zeros([round(steps), N, m]);
    p_history = zeros([round(steps), N, m]);
    qr_history = zeros([round(steps), N, m]);
    pr_history = zeros([round(steps), N, m]);
    de_history = zeros([round(steps), N, m]);
    
    %system evolution
    for t=0:tcyc:tmax-tcyc
        t_i=round(t/tcyc)+1;
        if t_i==1 || mod(t_i+1, u_steps / updates) == 0
            %we will call the below plotting loop every some odd iterations of the
            %system evolution to see how it's doing
            if m == 1 || m == 2 || m == 3
                figure;
                hold on;
            end

            if m == 1
                subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),zeros(N,1));
                subplot(subplotRows,subplotCols,subplotCount), scatter(qd(:,1)+pd(:,1).*t,zeros(N,1),'r');
                scatter(qd(:,1)+pd(:,1).*t,zeros(N,1),r,'g');
                scatter(qd(:,1)+pd(:,1).*t,zeros(N,1),d,'k');
            elseif m == 2
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
            for i=1:N
                neighbors_i = neighborsSpatialLattice(i, q(i,:), q, r, d, delta); %exclude self here if desired q_js: cat(1,q(1:i,:),q(i+1:size(q,1),:))
                %spatial_neighbors(:,:,i) = neighbors_i;

                %TODO: speed this up by removing the for loop and doing
                %everything via index operations after generating the neighbors
                for j=1:size(neighbors_i)
                    %prune data
                    %if j > i
                    %    continue;
                    %end
                    
                    if m == 1
                        plot([q(i,1) q(neighbors_i(j),1)], [0 0]);
                    elseif m == 2
                        plot([q(i,1) q(neighbors_i(j),1)], [q(i,2) q(neighbors_i(j),2)]);
                    elseif m == 3
                        plot3([q(i,1) q(neighbors_i(j),1)], [q(i,2) q(neighbors_i(j),2)], [q(i,3) q(neighbors_i(j),3)]);
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
        uGradient_history(t_i,:,:) = uGradient(:,:);
        uConsensus_history(t_i,:,:) = uConsensus(:,:);
        uGamma_history(t_i,:,:) = uGamma(:,:);
        
        %reinitialze all controls (not dependent upon past control value)
        u = zeros(N, m);
        uGradient = zeros(N, m);
        uConsensus = zeros(N, m);
        uGamma = zeros(N, m);
        
        %compute control (based on state vector, possibly delayed, etc)
        for i=1:N
            js = neighborsSpatial(i, q(i,:), q, r, d);
            %js = neighborsSpatialLattice(i, q(i,:), q, r, d, delta);
%            if size(js,1) > -1
                for j=1:size(js,1)
                    %TODO: verify equations and all functions
                    %u(i,:) = u(i,:) + phi_a(sig_norm ( q(j,:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * (((q(j,:) - q(i,:)) ) / (sqrt(1 + epsilon * ((norm(q(j,:) - q(i,:), 2))^2))));
                    %u(i,:) = u(i,:) + (a_ij(q(i,:), q(j,:), r_sig, ha,epsilon) * ((p(j,:) - p(i,:))));
                    uGradient(i,:) = uGradient(i,:) + phi_a(sig_norm ( q(js(j),:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * n_ij(q(i,:), q(js(j),:), epsilon);
                    uConsensus(i,:) = uConsensus(i,:) + (a_ij(q(i,:), q(js(j),:), r_sig, ha, epsilon) * ((p(js(j),:) - p(i,:))));
                end
%             else
%                 %no neighbors yet: have to compute something
%                 %TODO: fix, not right (always 0 obviously)
%                 u(i,:) = u(i,:) + phi_a(sig_norm ( q(i,:) - q(i,:), epsilon ), r_sig, d_sig, ha, a, b) * (((q(i,:) - q(i,:)) ) / (sqrt(1 + epsilon * norm(q(i,:) - q(i,:), 2))^2));
%                 u(i,:) = u(i,:) + (a_ij(q(i,:), q(i,:), r_sig, ha, epsilon) * ((p(i,:) - p(i,:))));
%             end
            
            %add gamma goal term
            %u(i,:) = u(i,:) - c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));
            uGamma(i,:) = -c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));
        end
        
        %sum all forces
        %uGamma(:,:) = zeros(N,m)
        u(:,:) = uGradient(:,:) + uConsensus(:,:) + uGamma(:,:);

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
            de_history(t_j,:,:) = de(:,:);
            
            %TODO: verify correctness of this as a solution
            %      in general, can we use this solution if the control
            %      is nonlinear but numerical?
            q = q + p.*(tcyc/tdiv) + u.*((tcyc/tdiv)^2);
            p = p + u.*(tcyc/tdiv);
            de = deviationEnergy(q,r,d,delta);
            
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
            
            if m == 1
                plot(time_ctrl,u_history(:,i,1),'b--');
                plot(time_ctrl,uGradient_history(:,i,1),'r--');
                plot(time_ctrl,uConsensus_history(:,i,1),'k--');
                plot(time_ctrl,uGamma_history(:,i,1),'g--');
                plot(time_traj,q_history(:,i,1) - qr_history(i,1),'c-.');
                plot(time_traj,p_history(:,i,1) - pr_history(i,1),'m-.');
                %plot(time_traj,de_history(:,i,1),'c');
                legend('u', 'uGradient', 'uConsensus', 'uGamma', 'q-qr', 'p-pr');
                %legend('u', 'uGradient', 'uConsensus', 'uGamma', 'de');
            elseif m == 2
                %plot(time_ctrl,u_history(:,i,1),'b--');
                %plot(time_ctrl,u_history(:,i,2),'c--');
                plot(time_ctrl,sqrt(u_history(:,i,1).^2 + u_history(:,i,2).^2),'b--');
                plot(time_ctrl,sqrt(uGradient_history(:,i,1).^2 + uGradient_history(:,i,2).^2),'r--');
                plot(time_ctrl,sqrt(uConsensus_history(:,i,1).^2 + uConsensus_history(:,i,2).^2),'k--');
                plot(time_ctrl,sqrt(uGamma_history(:,i,1).^2 + uGamma_history(:,i,2).^2),'g--');
                %plot(time_traj,q_history(:,i,1) - qd(i,1),'r:');
                %plot(time_traj,q_history(:,i,2) - qd(i,1),'m:');
                %plot(time_traj,p_history(:,i,1) - pd(i,1),'k:');
                %plot(time_traj,p_history(:,i,2) - pd(i,1),'b:');
                legend('u', 'uGradient', 'uConsensus', 'uGamma');
            elseif m == 3
                %use only norms here or we'll get too busy
                plot(time_ctrl,sqrt(u_history(:,i,1).^2 + u_history(:,i,2).^2 + u_history(:,i,3).^2),'b--');
                plot(time_ctrl,sqrt(uGradient_history(:,i,1).^2 + uGradient_history(:,i,2).^2 + uGradient_history(:,i,3).^2),'r--');
                plot(time_ctrl,sqrt(uConsensus_history(:,i,1).^2 + uConsensus_history(:,i,2).^2 + uConsensus_history(:,i,3).^2),'k--');
                plot(time_ctrl,sqrt(uGamma_history(:,i,1).^2 + uGamma_history(:,i,2).^2 + uGamma_history(:,i,3).^2),'g--');
                legend('u', 'uGradient', 'uConsensus', 'uGamma');
            end
        end
    end

    %average position and velocity for moving reference frame
    q
    p
    qd
    pd
    de
    qc=Ave(q)
    pc=Ave(p)
end


%1: fix lattice: goals independently in grid with multiple objectives

%2: 2 particle set up with one goal: how are these aligned? simplify
%functions enough for proofs

%3: obstacles

%crash failure: not moving, just obstacle

%stuck at failure: moving and growing in size: nodes know where it is and
%how fast it is going

%    => what is periodicity with which a node has to be able to update its
%       neighbors in order to be able to do mitigation in time

%   what if neighbor fails at max acc? can never go over it: how to handle
%   this?  when to overtake node and when to just follow it?

%4: failure with stuck value, starting in lattice

%5: one node slowly updating, everyone else ignores eventually and treats
%it as an obstacle: growing larger and larger obstacle as uncertainty is
%compounded over and over again

%john sicyclis: most general necessary and sufficient: on asynchronous
%iterated systems

%mani chandy's paper was like this going towards message passing systems




%when have we reached the goal?
%node is said to reach a goal if it's \alpha-lattice has reached the goal
%  may have to fine tune this definition



%1.a obstacles
%first: one obstacle

%1.b delayed updates?

%detection of failures: assume taking care of (node knows when it is going
%   over obstacle)
%avoidance: mitigation of failures

%2 piecewise psi function with 3-segments: 


%3 safety (minimum gap) and liveness/progress (convergence to goal): we do
%want to prove this


%4 work with their convergence (hard) and prove safety from there?
%problems arising from their synchronous assumptions?

%5 simplify everything and prove safety (easy) and convergence (hard)


%their algorithm may work for our model, but our model is inherently
%different

%how to define our energy function? tuple of distance-to-goal and deviation
%energy from \alpha-lattice?

%need velocity in goal term? oscillations? do oscillations provide a
%problem?

%relax condition such that nodes reach goal, maybe in formation, maybe not