function [ out ] = flocking(N, m, coord_min, coord_max, r_comm, r_lattice, r_safety, tdiv, tmax, updates, plotControls, tswitch)
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
    c1gamma=0.25;
    c2gamma=0.1;
    c1beta=0.25;
    c2beta=0.1;
    
    epsilon = 0.1;
    a=5;
    b=5;
    ha=0.2;
    hb=0.9;

    kappa = r_comm / r_lattice;

    v_max = 5;  %maximum velocity

    a_max = 100;    %maximum acceleration

    Tc=(min(2*v_max, (r_comm / (2*v_max / a_max))))/a_max %need to get r_comm into the calculation for Tc
    %Tc=0.01;
    %Tc=sqrt(r_safety/v_max*(v_max/a_max));
    %Tc=0.001;
    fc=1/Tc;
    t_a = 2*v_max / a_max
    t_v = r_safety / (2*v_max) 
    %r_init = r_safety + Tc*2*v_max + (Tc^2)*a_max + Tc*(2*v_max/a_max) + Tc*(r_safety/(2*v_max))
    %r_init = r_safety/(2*v_max/a_max)
    r_init = r_safety + max(2*v_max*((t_a)), v_max * t_v)
    r_lattice=r_init;
    
    delta=r_safety/5;
    
    r_lattice_prime = 0.6 * r_lattice;
    r_comm_prime = 1.2 * r_lattice_prime;
    
    r_comm_sig    = sig_norm(r_comm, epsilon);
    r_lattice_sig = sig_norm(r_lattice, epsilon);
    r_safety_sig  = sig_norm(r_safety, epsilon);
    
    delay_min = - (Tc / 4);
    delay_max = (Tc / 4);

    %generate velocity matrix
    %p = ones(N, m); %start from 1 velocity
    p = zeros(N, m)
    
    %p = v_max + (-v_max - v_max).*rand(N, m)
    
    %generate state matrix randomly such that no nodes are already
    %neighbors
%     q=zeros(N, m);
%     for i=1:1:N
%         q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);
%         
%         for j=1:1:i %only need to go up to i's already initialized
%             if i ~= j
%                 %ensure no vertices start as neighbors, that is 
%                 %(norm(q_i,q_j,2) < r) != true
%                 while size(neighborsSpatial(i, q(i,:), q, r_init, r_init),1) > 0
%                     q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);
%                 end
%             end
%         end
%     end
    
    %q=[0:(r_lattice-delta):(r_lattice-delta)*(N-1)]';
    
    %q = [5; 50; 90; 120; 150; 180; 210; 230; 250; 260]
    q = [1:N]'*r_init
    %q = [5; 5+r_safety+((r_init-r_safety)/2); 19; 26; 40; 47; 55; 62; 70; 77]
    %p = [v_max; -v_max; v_max; v_max; v_max; -v_max; v_max; -v_max; -v_max; v_max]
    %p = [v_max; -v_max; v_max; v_max; v_max; -v_max; v_max; -v_max; -v_max; v_max;
    %     v_max; -v_max; v_max; v_max; v_max; v_max; v_max: v_max: v_max; -v_max; -v_max; -v_max]
    %p = [v_max; -v_max];
    %v_max; v_max; v_max; v_max; v_max; v_max; -v_max; v_max; 
        %-v_max; v_max; v_max; v_max; v_max; v_max; v_max; v_max; v_max; -v_max]
    
    %q = [5; 5+r_init]
    %p = [v_max; -v_max]
    
    for c0=1:N
        for c1=1:m
            if (rand(1) > 0.5)
                p(c0,c1) = v_max;
            else
                p(c0,c1) = -v_max;
            end
        end
    end
    p
    
    
    
    size(q)
    size(p)
%    return

    spatial_neighbors = zeros(N, N, N);

    tcyc=Tc;
    time_ctrl=[0:tcyc:tmax-tcyc]';       %vector of control cycle times
    time_traj=[0:tcyc/tdiv:tmax]';  %vector of trajectories
    u_steps=round(((tmax-tcyc)/tcyc))+1;   %number of control update steps
    steps=tmax/(tcyc/tdiv);         %total numbers of trajectory steps

    %goals for gamma control term
    qd=zeros(N, m);         %preallocate
    %pd=ones(N, m);
    pd=zeros(N, m); %end with 0 velocity
    
    %initial control value
    %u = (vel_max + (vel_min - vel_max).*rand(N, m));
    u = zeros(N, m);
    uGradient = zeros(N, m);
    uConsensus = zeros(N, m);
    uGamma = zeros(N, m);
    uNew = zeros(N, m);

    %generate an equally spaced (by r_lattice) rectangular grid goal
    %TODO: generalize for m-dimensional (use cat and m-th root of N [need m loops?])
    if m == 1
        for i=1:N
            qd(i,:) = ((i-1)*r_lattice) + coord_max;
        end
        %scatter(qd(:,1),zeros(N,1),[1:N]');
    elseif m == 2
        for i=1:ceil(sqrt(N))
            for j=1:ceil(sqrt(N))
                qd((i-1)*ceil(sqrt(N))+j,:) = [((i-1)*r_lattice) ((j-1)*r_lattice)] + [coord_max coord_max];
            end
        end
        %scatter(qd(:,1),qd(:,2));
    elseif m == 3
        for i=1:ceil(N^(1/3))
            for j=1:ceil(N^(1/3))
                for k=1:ceil(N^(1/3))
                    qd((i-1)*ceil(N^(2/3))+(j-1)*ceil(N^(1/3))+k,:) = [((i-1)*r_lattice) ((j-1)*r_lattice) ((k-1)*r_lattice)] + [coord_max coord_max coord_max];
                end
            end
        end
        %scatter(qd(:,1),qd(:,2),qd(:,3));
    end
    qd=qd(1:N,:); %shrink to maximum N

    qd=ones(N,m).*coord_max*(2+N);
    pd=zeros(N,m);
    
    qr=qd;
    pr=pd;
    
    %start at goal
    %q=qd;
    %p=pd;
    
    %obstacles: k obstacles with radius Rk, centered at yk
    %Ms = [150 150; 30 100; 25 25]; %from paper, last row is radius, first rows are x,y,z,... positions
    %Ms = [100 110 120 130 150 160; 20 60 40 -20 40 0; 10 4 2 5 5 3]
    %R_k = Ms(size(Ms,1),:);
    %for i=1:size(Ms,2)
    %    for j=1:size(Ms,1)-1
    %        y_k(i,j) = Ms(j,i);
    %    end
    %end
    
    %de = deviationEnergy(q,r_comm,r_lattice,delta)
    %pne = proximityNetEdges(q,r_comm,r_lattice)
    
    subplotRows=1;%ceil(sqrt(updates));
    subplotCols=1;%floor(sqrt(updates));
    subplotCount=1;

    %control and state variables saved at points in time over evolution
    u_history = zeros([round(u_steps), N, m]);
    uGradient_history = zeros([round(u_steps), N, m]);
    uConsensus_history = zeros([round(u_steps), N, m]);
    uGamma_history = zeros([round(u_steps), N, m]);
    uNew_history = zeros([round(u_steps), N, m]);
    q_history = zeros([round(steps), N, m]);
    p_history = zeros([round(steps), N, m]);
    qr_history = zeros([round(steps), N, m]);
    pr_history = zeros([round(steps), N, m]);
    %de_history = zeros([round(steps), N, m]);
    
    %uPeriod = (1:N)'.*0.01 %different update periods for all particles
    uPeriod = ones(N,1)*Tc %same update period for all particles
    %uPeriod(1) = Tc*5; %make one node update its control slowly

    %add arbitrary delay to each control update
    %uOffset = delay_max + (delay_min - delay_max).*rand(N, 1)
    tadd = tcyc / tdiv;
    uOffset = zeros(N,1);
    for i=1:N
        to=round(rand(1,1)*tdiv);
        for j=1:to
            %uOffset(i) = uOffset(i) + tadd;
        end
    end
    uOffset

    switches = 0; %initialize to 0 (before any goal updates occur)

    %system evolution
    for t=0:tcyc:tmax-tcyc
        t_i=round(t/tcyc)+1;

        %update goals after some time
        %TODO: change to make new goals after reaching current goal (allows
        %      for arbitrary length system evolution that makes sense)
%         if t > tswitch && switches == 0
%             %set a new goal
%             qd = qd.^2 + ones(N, m)*100;
%             pd = pd.^2 + ones(N, m)*10;
%             qr = qd
%             pr = pd
%             switches = switches + 1;
%         end
%         
%         if t > (tswitch*3) && switches == 1
%             %set a new goal
%             qd = qd*2 + ones(N, m)*10;
%             pd = pd*2 + ones(N, m);
%             qr = qd
%             pr = pd
%             switches = switches + 1;
%         end
        if t_i==1 || mod(t_i+1, u_steps / updates) == 0
            %we will call the below plotting loop every some odd iterations of the
            %system evolution to see how it's doing
            if m == 1 || m == 2 || m == 3
                figure;
                hold on;
            end

            if m == 1
%                 isSafety = zeros(N,1);
%                 isLattice = zeros(N,1);
%                 isComm = zeros(N,1);
%                 isElse = zeros(N,1);
%                 qSafety = zeros(N,1);
%                 qLattice = zeros(N,1);
%                 qComm = zeros(N,1);
%                 qElse = zeros(N,1);
% 
%                 %locate nodes satisfying different radius constraints
%                 %(communication, lattice, and safety) and plot in various
%                 %colors so we can more easily see system evolution
%                 %TODO: would be nice to have a bubble plot with radius of
%                 %      circles the size of the safety region (and also
%                 %      additionally interaction region if different
%                 %      radius), but Matlab's bubble plotting seems rather
%                 %      limited as the radius appears to be a relative size,
%                 %      not an actual point value
%                 iSafety=0;
%                 iLattice=0;
%                 iComm=0;
%                 iElse=0;
%                 for i=1:N
%                     done = 0;
%                     bad = 0;
%                     %search over already bad nodes and don't add again
%                     for j=1:N
%                         if i == isSafety(j) || i == isLattice(j) || i == isComm(j) || i == isElse(j)
%                             done = 1;
%                             break;
%                         end
%                     end
%                     
%                     if done == 1
%                         continue; %already added this as bad
%                     end
%                     
%                     for j=1:N
%                         if i == isSafety(j) && (i ~= j)
%                             bad = 1;
%                             break;
%                         elseif (norm(q(i,1) - q(j,1), 2) <= r_safety) && (i ~= j)
%                             %this needs to occur before other as anything
%                             %that satisfies safety would satisfy lattice
%                             %and comm
%                             qSafety(iSafety + 1) = q(i,:);
%                             isSafety(iSafety + 1) = i; %save indexes
%                             iSafety = iSafety + 1;
%                             bad = 1;
%                         elseif (norm(q(i,1) - q(j,1) - r_lattice, 2) <= delta) && (i ~= j)
%                             %TODO: change <= r_lattice to |qi - qj - r_lattice| <= delta
%                             %      to allow for quasi-lattices
%                             qLattice(iLattice + 1) = q(i,:);
%                             isLattice(iLattice + 1) = i; %save indexes
%                             iLattice = iLattice + 1;
%                             bad = 1;
%                         elseif (norm(q(i,1) - q(j,1), 2) <= r_comm) && (i ~= j)
%                             qComm(iComm + 1) = q(i,:);
%                             isComm(iComm + 1) = i; %save indexes
%                             iComm = iComm + 1;
%                             bad = 1;
%                         end
%                     end
%                     
%                     if bad == 0
%                         %if wasn't bad, add to good
%                         qElse(iElse + 1) = q(i,:);
%                         iElse = iElse + 1;
%                     end
%                 end
%                 qElse = qElse(1:iElse); %shrink
%                 qComm = qComm(1:iComm); %shrink
%                 qLattice = qLattice(1:iLattice); %shrink
%                 qSafety = qSafety(1:iSafety); %shrink
%                 qElse
%                 qComm
%                 qLattice
%                 qSafety
%                 scatter(qElse(:,1),zeros(iElse,1),'k');
%                 scatter(qComm(:,1),zeros(iComm,1),r_comm,'b');
%                 scatter(qLattice(:,1),zeros(iLattice,1),r_lattice,'g');
%                 scatter(qSafety(:,1),zeros(iSafety,1),r_safety,'r');
%                 %subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),zeros(N,1),'b'); 
%                 scatter(qr(:,1),zeros(N,1),'g'); %plot goal
                subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),zeros(N,1),'b');
                scatter(qr(:,1),zeros(N,1),'g'); %plot goal
            elseif m == 2
                subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),q(:,2));
                subplot(subplotRows,subplotCols,subplotCount), scatter(qr(:,1),qr(:,2),'r');
                %out=[[qd(:,1) (qd(:,1)+pd(:,1))] [qd(:,2) (qd(:,2)+pd(:,2))]];
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplotCount = subplotCount + 1;
            elseif m == 3
                subplot(subplotRows,subplotCols,subplotCount), scatter3(q(:,1),q(:,2),q(:,3));
                subplot(subplotRows,subplotCols,subplotCount), scatter3(qr(:,1),qr(:,2),qr(:,3),'r');
                %subplotCount = subplotCount + 1;
            end

            %generate neighbors for all nodes and plot
            %TODO: have this simply create a data structure and then plot that
            %      What will that data structure look like?  Rather complicated
            %      multi-dimensional object without consistent number of elements
            %TODO: move to function
            for i=1:N
                neighbors_i = neighborsSpatialLattice(i, q(i,:), q, r_comm, r_lattice, delta); %exclude self here if desired q_js: cat(1,q(1:i,:),q(i+1:size(q,1),:))
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

        for t_j=(t_i-1)*tdiv+1 : 1 : tdiv+(t_i-1)*tdiv+1
            tt=t_j*(tcyc/tdiv);

            %compute control (based on state vector, possibly delayed, etc)
            for i=1:N
                if t_j == 1 || (mod(tt, (uPeriod(i) + uOffset(i))) == 0)
                    %tt
                    %u_i = u_i^\alpha + u_i^\gamma
                    %n_ij = ((q_j - q_i ) / (sqrt(1 + epsilon * norm(q_j - q_i, 2))^2));

                    %rho_h(z) = 1                                    if z \in [0,h)
                    %           1/2 * (1 + cos(pi * (z - h)/1 -h))   if z \in [h,1]
                    %           0                                    else
                    %for h \in (0,1)

                    %a_ij(q) = rho_h( sig_norm ( q_j - q_i ) / r_comm_sig )
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
                    uNew_history(t_i,:,:) = uNew(:,:);

                    %reinitialze this nodes controls (not dependent upon past control value)
                    u(i,:) = zeros(1,m);
                    uGradient(i,:) = zeros(1,m);
                    uConsensus(i,:) = zeros(1,m);
                    uBetaGradient(i,:) = zeros(1,m);
                    uBetaConsensus(i,:) = zeros(1,m);
                    uGamma(i,:) = zeros(1,m);
                    uNew(i,:) = zeros(1,m);

                    js = neighborsSpatial(i, q(i,:), q, r_comm, r_lattice);
                    betas = neighborsSpatial(i, q(i,:), q, r_comm_prime, r_lattice_prime);
                    %js = neighborsSpatialLattice(i, q(i,:), q, r_comm, r_lattice, delta);
        %            if size(js,1) > -1
                        for j=1:size(js,1)
                            %TODO: verify equations and all functions
                            %u(i,:) = u(i,:) + phi_a(sig_norm ( q(j,:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * (((q(j,:) - q(i,:)) ) / (sqrt(1 + epsilon * ((norm(q(j,:) - q(i,:), 2))^2))));
                            %u(i,:) = u(i,:) + (a_ij(q(i,:), q(j,:), r_comm_sig, ha,epsilon) * ((p(j,:) - p(i,:))));
                            %old: uGradient(i,:) = uGradient(i,:) + phi_a(sig_norm ( q(js(j),:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * n_ij(q(i,:), q(js(j),:), epsilon);
                            %old: uConsensus(i,:) = uConsensus(i,:) + (a_ij(q(i,:), q(js(j),:), r_comm_sig, ha, epsilon) * ((p(js(j),:) - p(i,:))));
                            uGradient(i,:) = uGradient(i,:) + phi_a(sig_norm ( q(js(j),:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * n_ij(q(i,:), q(js(j),:), epsilon);
                            uConsensus(i,:) = uConsensus(i,:) + (a_ij(q(i,:), q(js(j),:), r_comm_sig, ha, epsilon) * ((p(js(j),:) - p(i,:))));
                        end
                        
                        for j=1:size(betas,1)
                            uBetaGradient(i,:) = uBetaGradient(i,:) + phi_a(sig_norm ( q(betas(j),:) - q(i,:), epsilon ), r_comm_prime, r_lattice_prime, ha, a, b) * n_ij(q(i,:), q(betas(j),:), epsilon);
                            uBetaConsensus(i,:) = uBetaConsensus(i,:) + (a_ij(q(i,:), q(betas(j),:), r_comm_prime, ha, epsilon) * ((p(betas(j),:) - p(i,:))));
                        end
        %             else
        %                 %no neighbors yet: have to compute something
        %                 %TODO: fix, not right (always 0 obviously)
        %                 u(i,:) = u(i,:) + phi_a(sig_norm ( q(i,:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * (((q(i,:) - q(i,:)) ) / (sqrt(1 + epsilon * norm(q(i,:) - q(i,:), 2))^2));
        %                 u(i,:) = u(i,:) + (a_ij(q(i,:), q(i,:), r_comm_sig, ha, epsilon) * ((p(i,:) - p(i,:))));
        %             end

                    %add gamma goal term
                    %u(i,:) = u(i,:) - c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));
                    %old: uGamma(i,:) = -c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));

                    %uGamma(i,:) = -c1gamma*sigma_1((q(i,:) - qr(i,:))) - c2gamma*(p(i,:) - pr(i,:));
                    uGamma(i,:) = -c1gamma*((q(i,:) - qr(i,:))) - c2gamma*(p(i,:) - pr(i,:));

                    %sum all forces
                    %uGamma(:,:) = zeros(N,m)
                    %olfati-saber control
                    %u(i,:) = uGradient(i,:) + uConsensus(i,:) + uGamma(i,:) + c1beta*uBetaGradient(i,:) + c2beta*uBetaConsensus(i,:);
                    
                    %1: check all nodes to see if they see anyone on their
                    %   right, then have them go to goal with max
                    %   acceleration
                    %   for implementation, just check THIS node (i) and
                    %   see if this is the case
                    %2: find nearest right hand neightbor and apply 
                    %   u_i = a * (xhat_i - x_i - s) + b * (vhat_i - v_i)
                    %   where s is lattice separation (r_lattice)

                    uNew(i,:) = ((qd(i,:) - q(i,:)) / norm(qd(i,:) - q(i,:), 2)) * a_max;
                    dist_min_right = inf;
                    dist_min_left = inf;

                    %overwrite control if necessary
                    for j=1:N
                        if (norm(q(j) - q(i),2) <= r_comm) && ((q(j)) > (q(i))) && (j ~= i) && (norm(q(j) - q(i),2) <= dist_min_right)
                            %uNew(i,:) = aS * (q(j,:) - q(i,:) - r_lattice + delta) + bS * (p(j,:) - p(i,:));
                            uNew(i,:) = (1/Tc) * (q(j,:) - q(i,:) - r_init) + (1/Tc) * (p(j,:) - p(i,:));
                            %uNew(i,:) = a_max;
                            dist_min_right = norm(q(j) - q(i),2);
                            %if norm(q(j) - q(i),2) <= (r_lattice - delta)
                            if norm(q(j) - q(i),2) <= r_init
                                uNew(i,:) = -a_max;
                            end
%                         elseif (norm(q(j) - q(i),2) <= r_comm) && ((q(j)) < (q(i))) && (j ~= i) && (norm(q(j) - q(i),2) <= dist_min_left)
%                             %slow the lead node down to form flocking
%                             dist_min_left = norm(q(j) - q(i),2);
%                             uNew(i,:) = -a_max;
%                             %uNew(i,:) = - (1/Tc) * (q(j,:) - q(i,:) - r_init) - (1/Tc) * (p(j,:) - p(i,:));
%                             %uNew(i,:) = -(aS/4) * (q(j,:) - q(i,:) - r_init) - bS * (p(j,:) - p(i,:));
%                             if norm(q(j) - q(i),2) <= r_init
%                                 uNew(i,:) = a_max;
%                             end
                        %elseif (norm(q(j) - q(i),2) <= r_comm) && ((q(j)) < (q(i))) && (j ~= i) && (norm(q(j) - q(i),2) <= dist_min)
                        %    uNew(i,:) = -aS * (q(j,:) - q(i,:) - r_lattice) - bS * (p(j,:) - p(i,:));
                        %    dist_min = norm(q(j) - q(i),2);
                        end
                        
                        %start from conditions: 0 velocity, close to one
                        %another (within comm distance and almost safety,
                        %such that radius to form lattice is between these
                        %two)
                        
                        %add condition such that if neighbor to right is
                        %too close, wait instead of going in opposite
                        %direction
                    end
                    u(i,:) = uNew(i,:);
                end
            end
            
            %saturate the controls (no effect until out of neighbors loop
            %above, so might as well be efficient and do it once for
            %everything)
            u = sign(u).*(min(abs(u),a_max));
            
            %run system evolution
            
            %saturate the velocity
            %sign(p)
            %(min(abs(p),v_max))
            p = sign(p).*(min(abs(p),v_max));
            
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
            %de_history(t_j,:,:) = de(:,:);
            
            %TODO: verify correctness of this as a solution
            %      in general, can we use this solution if the control
            %      is nonlinear but numerical?
            q = q + p.*(tcyc/tdiv) + u.*((tcyc/tdiv)^2);
            p = p + u.*(tcyc/tdiv);
            %de = deviationEnergy(q,r_comm,r_lattice,delta);
            
            %gamma agent
            qr=qr + pr.*(tcyc/tdiv) + fr(qr, pr).*((tcyc/tdiv)^2);
            pr=pr + fr(qr, pr).*(tcyc/tdiv);
            
            %beta agent
            %mu = R_k / norm(q_i - y_k, 2);
            %a_k = (q_i - y_k) / norm(q_i - y_k, 2);
            %P = I - a_k * a_k';
            %qhat_ik = mu * q_i + (1 - mu) * y_k;
            %phat_ik  = mu * P * p_i;
        end
    end
    
    out = u_history;

    %plot controls over time
    if plotControls >= 1
        if plotControls == 3
            %plot distance of each node so we can see how their lattice
            %evolves
            figure;
            hold on;
            for i=1:N
                if mod(i,2)==0
                    plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'r');
                else
                    plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'b');
                end
            end
            legend('q1');

            if m >= 2
                figure;
                hold on;
                for i=1:N
                    plot(time_traj,q_history(:,i,2),'m-.');
                    legend('q2');
                end
            end
        end

        if plotControls == 2
            for i=1:N
                figure;
                hold on;

                if m == 1
                    plot(time_ctrl,u_history(:,i,1),'b--');
                    %plot(time_ctrl,uGradient_history(:,i,1),'r--');
                    %plot(time_ctrl,uConsensus_history(:,i,1),'k--');
                    %plot(time_ctrl,uGamma_history(:,i,1),'g--');
                    plot(time_ctrl,uNew_history(:,i,1),'m--');
                    plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'c-.');
                    plot(time_traj(1:size(p_history(:,i,1))),p_history(:,i,1),'g-.');
                    %plot(time_traj,q_history(:,i,1) - qr_history(i,1),'c-.');
                    %plot(time_traj,p_history(:,i,1) - pr_history(i,1),'m-.');
                    %plot(time_traj,de_history(:,i,1),'c');
                    %legend('u', 'uGradient', 'uConsensus', 'uGamma', 'q-qr', 'p-pr');
                    %legend('u', 'uGradient', 'uConsensus', 'uGamma', 'de');
                    legend('u', 'uNew', 'q', 'p');
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

                    if plotControls == 2
                        figure;
                        hold on;
                        plot3(time_ctrl,u_history(:,i,1),u_history(:,i,2),'b--');
                        plot3(time_ctrl,uGradient_history(:,i,1),uGradient_history(:,i,2),'r--');
                        plot3(time_ctrl,uConsensus_history(:,i,1),uConsensus_history(:,i,2),'k--');
                        plot3(time_ctrl,uGamma_history(:,i,1),uGamma_history(:,i,2),'g--');
                        plot3(time_traj,q_history(:,i,1) - qd(i,1),q_history(:,i,2) - qd(i,2),'r:');
                        plot3(time_traj,p_history(:,i,1),p_history(:,i,2),'k:');
                        legend('u', 'uGradient', 'uConsensus', 'uGamma', 'q - qd', 'p');
                    end
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
    end

    %average position and velocity for moving reference frame
    q
    p
    qd
    pd
    %de
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








%split rejoin setup
%q0: [-40, 80]^2
%p0: [0, 0]
%group objective, static gamma agent: 
%qd=[200,30]'
%pd=[5,0]'
%c1alpha < c1gamma < c1beta
%c2nu = 2 * sqrt(c1nu)  (for all beta, gamma, alpha constants)
%
%obstacles:

%Ms = [100 110 120 130 150 160; 20 60 40 -20 40 0; 10 4 2 5 5 3]

%squeezing manaeuver:
%n=150
%q0: [0, 120]^2
%p0: [0, 0]^2
%group objective: static gamma agent
%qd=[230, 60]'
%pd=[6, 0]'
%c1alpha < c1gamma < c1beta
%c2nu = 2 * sqrt(c1nu) (for all beta, gamma, alpha constants)

%obstacles:
%Ms = [150 150; 30 100; 25 25]


%Ms is defined as:
%The set of l spherical obstacles is specified as the
%(m + 1) by l matrix Ms where each column of Ms is the
%vector col(y_k, R_k) \in R^(m+1)

%radius Rk center at yk
%mu = R_k \ norm(q_i - y_k, 2)
%a_k = (q_i - y_k) \ norm(q_i - y_k, 2)
%P = I - a_k * a_k'
%qhat_ik = mu * q_i + (1 - mu) * y_k
%phat_ik  = mu * P * p_i