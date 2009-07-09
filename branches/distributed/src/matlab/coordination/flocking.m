function [ out ] = flocking(framework, N, m, coord_min, coord_max, r_comm, r_lattice, r_safety, tdiv, tmax, updates, plotControls, tswitch, delay, constrain)
%Flocking and Consensus Problem Simulations
%
%Taylor Johnson
%johnso99@NOSPAM-illinois.edu
%University of Illinois at Urbana-Champaign
%Coordinated Science Laboratory
%Department of Electrical and Computer Engineering
%
%       framework:  selects which model (controls, discrete-time vs.
%                   continuous time, etc.) to use
%                   0: model from Olfati-Saber, 2005
%                   1: new model based on Olfati-Saber
%                   2: model from Gazi and Passino, 2005
%       N:          number of agents/nodes/flock-members (>=1)
%       m:          number of dimensions (>=1)
%       coord_min:  will be removed, set to 0
%       coord_max:  will be removed, set to 1000 for example
%       r_comm:     communications radius (try 8.4)
%       r_lattice:  flocking/lattice/comfortable distance (try 7.0)
%       r_safety:   minimum distance to maintain safety (try 7.0)
%       tdiv:       number of times to allow the system to evolve between 
%                   control cycle updates (>= 1; for example, =1 is one
%                   division per control cycle is the same as only looking 
%                   at control cycle)
%       tmax:       time to run the system for (t_final)
%                   frameworks 0, 1: try ~20
%                   framework 2: try >= 1000
%       updates:    number of times to update the plots (>=1)
%       plotControls:   if 1, plots all controls as functions of time
%       tswitch:    used for multiple goals
%       delay:      enable asynchronism? 0 synchronous, 1 asynchronous
%       constrain:  enable control (and state) saturation constraints? 0 
%                   allows infinite valued controls and states, 1 
%                   constrains/saturates controls and states

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
    c1gamma=0.1;
    c2gamma=0.1;
    c1beta=0.25;
    c2beta=0.1;
    
    epsilon = 0.1;
    a=5;
    b=5;
    ha=0.2;
    hb=0.9;

    kappa = r_comm / r_lattice;
    
    if framework == 0
        system_type = 0; %0 continuous, 1 discrete-event system
        numLyap = 0;
        leadNode = 1;
        
        Tc = 0.01;
        r_init = r_comm;

        v_max = 100;     %maximum velocity
        a_max = 1000;     %maximum acceleration
        %a_max = inf;
        delta=r_lattice/5
    elseif framework == 1
        system_type = 0; %0 continuous, 1 discrete-event system
        numLyap = 0;
        leadNode = 1;

        v_max = 5;      %maximum velocity
        a_max = 10;     %maximum acceleration
        
        %invalid initial conditions
        if (r_comm <= r_safety)
            return;
        end

        t_a = 2*v_max / a_max %time to force v from v_max to -v_max (and vice-versa)
        t_v = r_safety / (2*v_max) %time to cover safety distance (doesn't make sense)
        %t_v = norm(r_comm - r_safety,2) / (2*v_max);
        %t_r = norm(r_comm - r_safety,2) / (2*v_max) %time to cover the difference between communications radius and safety radius

        %Tc = min([t_a, t_v, t_r])
        Tc = min([t_a, t_v])
        Tc = 0.01;
        r_a = 2*v_max * t_a %distance covered while fixing velocity
        r_v = 2*v_max * t_v %distance covered
        %r_r = 2*v_max * t_r
        r_d = 2*v_max * Tc %max distance covered over delay
        %r_init = r_safety + max([r_a, r_v, r_r]) + r_d
        %r_init = r_safety + max([r_a, r_v]) + r_d
        r_init = r_safety + r_a + r_d
        r_comm = r_init + r_d
        %r_init must be: r_safety + (max distance traveled before correcting
        %                control can stop [diff eq]) + (distance traveled
        %                during communications delay)
        r_lattice=r_init + r_d;
        delta=r_lattice/5
    elseif framework == 2
        system_type = 1; %0 continuous, 1 discrete-event system
        
        Tc = 1;
        
        alp_min = 0.001;
        alp_max = 0.999;
        eps = 0.25;
        dist = r_lattice;
        delta = r_lattice/5;
        
        jumpError = delta/2; %max error between any node before node 1 jumps
        
        step_max = 100000; %note, if we constrain the states like this, there will arrise a 
                           %stronger condition on the initial states: there must be a lower
                           %(safety) and upper (communications? kind of) bound on the 
                           %starting positions
        
        v_max = 1000;
        a_max = 1000;
        
        r_init = 5;
        
        flocking_weak_count = 0;
        
        kjump = 1;
        
        leadNode = ceil(N/3);
        
        numLyap = 3;
    elseif framework == 3 %inverted pendulum dynamics
        system_type = 0;
        leadNode = 1;
        delta=r_lattice/5;
        numLyap = 1;
        
        m = 4; %requires m = 4

        theta_min_deg=-30;  %degrees
        theta_max_deg=30;   %degrees
        theta_min_rad=deg2rad(theta_min_deg);   %radians
        theta_max_rad=deg2rad(theta_max_deg);   %radians
        %x_min=-0.7;         %meters
        %x_max=0.7;          %meters
        x_min=-0.2;         %meters
        x_max=0.2;          %meters
        xdot_min=-1;        %meters/second
        xdot_max=1;         %meters/second
        va_min=-4.95;       %volts
        va_max=4.95;        %volts
        
        v_max = xdot_max;

        Ts=0.020;
        Tc=0.020;   %20ms control period
        fs=1/Ts;
        fc=1/Tc;

        syms x1 x2 x3 x4;     %state variables
        syms x01 x02 x03 x04; %initial conditions
        syms a22 a23 a24 a42 a43 a44; %system matrix
        syms a22w a23w a24w a42w a43w a44w; %system matrix
        syms b2 b4; %control matrix
        syms Va;

        %physical parameters
        g=-9.8; %gravity: -9.8 m/s
        m=1;    %mass
        Ki=0.9; %torque constant
        Kb=0.9; %back-emf constant
        Kg=1;   %gear ratio
        r=2; %
        Ra=0.95;    %armature resistance
        Bm=1.95;    %
        Btheta=0.95;    
        Bbar=((Kg*Bm)/(r^2))+(((Kg^2)*Ki*Kb)/((r^2)*Ra));
        Bl=((Kg*Ki)/(r*Ra));
        M=5;    %mass
        Jm=5;   %moment of inertia
        l=0.5;  %length
        Mbar=m+M+(Kg*Jm)/(r^2); %effective mass
        Dl=(4*Mbar)-3*m;

        X0s=[x01 x02 x03 x04]';
        X=[x1 x2 x3 x4];
        U=[Va];
        %As=[0 1 0 0; 0 -a22 -a23 a24; 0 0 0 1; 0 a42 a43 -a44];
        As=[0 0 1 0; ; 0 0 0 1; 0 -a22 -a23 a24; 0 a42 a43 -a44];
        Asw=[0 1 0 0; 0 -a22w -a23w a24w; 0 0 0 1; 0 a42w a43w -a44w];
        %Bs=[0; b2; 0; -b4];
        Bs=[0; 0; b2; -b4];
        %a22=0.1;  a23=0.01;    a24=0.01;  a42=0.01;  a43=0.5;  a44=0.1;
        %a22w=(4*Bbar)/(Dl);  a23w=(3*m*g)/(Dl);   a24w=(6*Btheta)/(l*Dl);
        %a42w=(6*Bbar)/(l*Dl);    a43w=(6*Mbar*g)/(l*Dl);  a44w=(12*Mbar*Btheta)/(m*(l^2)*Dl);
        a22=2.75; a23=10.95; a24=0.0043; a42=28.58; a43=24.92; a44=0.044; %real matrix

        %Xreal=[x theta xdot thetadot]'

        %b2=0.1;   b4=0.5;
        %b2=(4*Bl)/(Dl); b4=(6*Bl)/(l*Dl);
        b2=1.94; b4=4.44; %real matrix
        %x1=x, x2=x', x3=theta, x4=theta'
        %x01=0.1; x02=0.2; x03=0.5; x04=-0.5;
        %x01=0.05; x02=0.05; x03=0.05; x04=0.05; %interesting plot: edge of safety
        %area
        %x01=0.05; x02=0.05; x03=-0.05; x04=0.05;
        %x01=0.1; x02=0.15; x03=0.05; x04=0.05;
        %x01=0.05; x02=0.05; x03=-0.1; x04=-0.1; %edge of larger safety
        %x01=-0.15; x02=-0.21; x03=0.4; x04=1.25;
        x01=-0.19; x02=-0.2418; x03=0.3; x04=1.41; %extrema safety
        %x01=-1; x02=-0.2581; x03=0.13; x04=1.17; %extrema experimental

        A=eval(As);
        Aw=eval(Asw);
        B=eval(Bs);
        X0=eval(X0s);

        %A=[37.62, 58.22, 17.87, 11.61; 58.22, 313.16, 69.36, 56.09; 17.87, 69.36, 29.81, 14.81; 11.61, 56.09, 14.81, 12.04];
        %A=[1 -0.00051281 0.017961 -0.0000026781; 0 1.0056 0.0046419 0.020029; 0 -0.049519 0.80322 -0.00043546; 0 0.55967 0.44824 1.0048;]

        %B=[0.0003618; -0.00082708; 0.034913; -0.079879];

        Q=eye(4); %positive definite Q
        syms p11 p12 p13 p14 p22 p23 p24 p33 p34 p44;
        P=[p11 p12 p13 p14; p12 p22 p23 p24; p13 p23 p33 p34; p14 p24 p34 p44];

        %linear state feedback control: V_a=KX=inv(-R)B'SX
        %A'S+SA-SBinv(R)B'S+D=0
        %R1=0.01
        %R2=0.01
        %D=eye(4)

        syms KS1 KS2 KS3 KS4 KB1 KB2 KB3 KB4 KE1 KE2 KE3 KE4;
        %KS=[7.6,  13.54, 42.85,  8.25];     %safety
        %KS=[6.0, 20.0, 60.0, 16.0];         %safety (real system)
        KS=[7.6,  42.85, 13.54,  8.25];     %safety
        KSs=[KS1, KS2, KS3, KS4];           %safety symbolic
        %KB=[3.16, 19.85, 69.92,  14.38];    %baseline
        %KB=[8.0, 32.0, 120.0, 12.0];        %baseline (real system)
        KB=[3.16, 69.92, 19.85,  14.38];    %baseline
        KBs=[KB1, KB2, KB3, KB4];
        %KE=[10.0, 27.72, 103.36, 23.04];    %experimental
        %KE=[5.7807, 42.2087, 14.0953, 8.6016];
        KE=[10.0, 103.36, 27.72, 23.04];    %experimental
        %KE=[10.0, 36.0, 140.0, 14.0];       %experimental (real system)
        KEs=[KE1, KE2, KE3, KE4];

        %from lmi
        Plmi=[38.3367   26.2506    8.6703    4.7515; 26.2506   56.5561   15.2589    8.6971; 8.6703   15.2589    6.6481    3.2088; 4.7515    8.6971    3.2088    1.8244];
        %Klmi=[7.2671   30.1022   11.7542    5.9111];
        %KS=Klmi;
        det(Plmi); %must be greater than 0 (shows P is positive definite)

        %Xdot=Abar*X
        AbarS=A+B*KS;
        AbarSs=As+Bs*KSs;
        AbarB=A+B*KB;
        AbarE=A+B*KE;

        det(AbarS'*Plmi+Plmi*AbarS); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)
        det(AbarB'*Plmi+Plmi*AbarB); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)
        det(AbarE'*Plmi+Plmi*AbarE); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)

        C=[0 1 0 1];
        sysS=ss(A,B*KS,C,0);
        sysB=ss(A,B*KB,C,0);
        sysE=ss(A,B*KE,C,0);
        %figure;
        %nyquist(sysS);
        %hold on;
        %nyquist(sysB);
        %nyquist(sysE);
        %legend('Safety', 'Baseline', 'Experimental');

        sysdS=ss(A,B*KS,C,0,Tc);
        sysdB=ss(A,B*KB,C,0,Tc);
        sysdE=ss(A,B*KE,C,0,Tc);
        %figure;
        %nyquist(sysdS);
        %hold on;
        %nyquist(sysdB);
        %nyquist(sysdE);
        %legend('Safety', 'Baseline', 'Experimental');

        %syms t;
        %Abarexpm=expm(AbarS*t);

        eigAbarS=eig(AbarS)
        eigAbarB=eig(AbarB)
        eigAbarE=eig(AbarE)

        %Lyapunov and Discrete Lyapunov Equations for each System
        P3=lyap(AbarS',Q)
        PD3=dlyap(AbarS',Q)
        P2=lyap(AbarB',Q)
        PD2=dlyap(AbarB',Q)
        P1=lyap(AbarE',Q)
        PD1=dlyap(AbarE',Q)

        SPtBad=0;
        SPdtBad=0;

        %set up for previous state (for switching)
        %let's run the entire system's 3 controllers concurrently, so that we can easily switch between them
        sc = 1; %safety
        bc = 2; %baseline
        ec = 3; %experimental
        swc = 4; %switched (properly, state-dependent)
        state(sc, 1:4) = X0;
        state(bc, 1:4) = X0;
        state(ec, 1:4) = X0;
        state(swc, 1:4) = X0;
        
        control(sc, 1) = 0;
        control(bc, 1) = 0;
        control(ec, 1) = 0;
        control(swc, 1) = 0;

        %Timing is crucial to our simulation here.  We use the given control cycle
        %time of Tc=20ms to control the main operation (that is, the control is 
        %updated every 20ms), then tdiv is the number of
        %subsamples to also model within this main period, so that we have smoother
        %trajectories within this larger control cycle.
        tcyc=Tc;

        %these are constants
        F=expm(A*Tc);
        Ft=expm(A*(tcyc/tdiv));
        PF=lyap(F,-Q);
        %F=A*Tc;
        syms tau;
        %G=int(expm(A*tau),tau,0,Tc);
        %G=A*(expm(A*Tc)-eye(4)); %int(expm(A*tau),tau,0,Tc)=A*expm(A*tau)|0,Tc=A*(expm(A*T)-expm(A*0))=A*(expm(A*T)-I)
        %[v,d]=eig(A);
        %G=eval(int((v*diag(exp(diag(d)))/v)*tau,0,Tc));
        G=pinv(A)*(expm(A*Tc)-eye(4));
        Gt=pinv(A)*(expm(A*(tcyc/tdiv))-eye(4));
        %G=inv(A)*((A*Tc)-eye(4));
        %from mathematica:
        %G=[0.02, 0.000199867, -1.32935*10^-8, 1.32535*10^-8; 4.38118*10^-19,
        %0.01998, -1.99204*10^-6, 1.98405*10^-6;7.78952*10^-21, 1.33201*10^-8,
        %0.0200007, 0.00019987; 1.47391*10^-20, 1.99737*10^-6, 0.0000999349, 0.0199807];

        cS=norm(P3*B,inf)/(min(eig(Q)));
        dS=norm(KS*A,inf)+norm(KS*B*KS,inf);
        tauS=1/(cS*dS);

        cB=norm(P2*B,inf)/(min(eig(Q)));
        dB=norm(KB*A,inf)+norm(KB*B*KB,inf);
        tauB=1/(cB*dB);

        cE=norm(P1*B,inf)/(min(eig(Q)));
        dE=norm(KE*A,inf)+norm(KE*B*KE,inf);
        tauE=1/(cE*dE);
    end
    
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
    if constrain ~= 0
        p = sign(p).*(min(abs(p),v_max));
    end
    
    q=zeros(N, m);
    
    if framework == 0
        %generate state matrix randomly such that no nodes are already neighbors
        for i=1:1:N
            q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);

            for j=1:1:i %only need to go up to i's already initialized
                if i ~= j
                    %ensure no vertices start as neighbors, that is 
                    %(norm(q_i,q_j,2) < r) != true
                    %while size(neighborsSpatial(i, q(i,:), q, r_init, r_init),1) > 0
                    %    q(i,:) = coord_max + (coord_min - coord_max).*rand(1, m);
                    %end
                end
            end
        end
    elseif framework == 1
        q(1:N/2) = [1:N/2]'*(r_init)
        q(N/2:N) = [N/2:N]'*(r_init+2) + 500
        
        for c0=1:N
            for c1=1:m
                if (rand(1) > 0.5)
                    p(c0,c1) = v_max; %worst case
                else
                    p(c0,c1) = -v_max;
                end
            end
        end
    elseif framework == 2
        for c0 = 1 : N
            for c1 = 1 : m
                if c0 == 1
                    q(c0,c1) = 0;
                else
                    q(c0,c1) = q(c0 - 1,c1) + eps + rand(1,1)*(dist*2.5 - eps);
                    %q(c0,c1) = q(c0 - 1,c1) + eps;
                    
                    %if mod(c0, 2) == 1
                    %    q(c0,c1) = q(c0 - 1,c1) + eps;
                    %else
                    %    q(c0,c1) = q(c0 - 1,c1) + dist*1.25;
                    %end
                    
                    %q(c0,c1) = q(c0 - 1) + (dist);
                    
                    if (abs(q(c0,c1) - q(c0 - 1,c1))) < eps
                        'bad safety initial condition'
                        q(c0,c1) = q(c0 - 1, c1) + eps; %needs to be above eps
                    end
                end
            end
        end
        q
        
        %q(1:N) = [1:N]'*(dist*3)

        %q(1:N/3) = [1:N/3]'*(r_init)
        %q(N/3:2*N/3) = [N/3:2*N/3]'*(r_init+2) + 25
        %q(2*N/3:N) = [2*N/3:N]'*(1) + 25
        
        for c0 = 1:m
            %q_goal(:,c0) = -25 + [0:N-1]'.*r_lattice; %based on node 1 as goal
            goal = -25;
            q_goal(:,c0) = goal + [-leadNode + 1: N - leadNode]'.*r_lattice;
            %q_goal(leadNode,c0) = goal;
            %q_goal(1:leadNode-1,c0) = goal + [0:leadNode-1]'.*r_lattice;
            %q_goal(leadNode+1:N,c0) = goal + [leadNode+1:N]'.*r_lattice;
        end
        leadNode
        N
        q_goal
        q_goal(leadNode)
    elseif framework == 3
        
        for c0 = 1 : N
            for c1 = 1 : m
                if c0 == 1
                    q(c0,c1) = 0;
                else
                    q(c0,c1) = q(c0 - 1,c1) + eps + rand(1,1)*(dist*2.5 - eps);
                    %q(c0,c1) = q(c0 - 1,c1) + eps;
                    
                    %if mod(c0, 2) == 1
                    %    q(c0,c1) = q(c0 - 1,c1) + eps;
                    %else
                    %    q(c0,c1) = q(c0 - 1,c1) + dist*1.25;
                    %end
                    
                    %q(c0,c1) = q(c0 - 1) + (dist);
                    
                    if (abs(q(c0,c1) - q(c0 - 1,c1))) < eps
                        'bad safety initial condition'
                        q(c0,c1) = q(c0 - 1, c1) + eps; %needs to be above eps
                    end
                end
            end
        end
        q
        
        %q(1:N) = [1:N]'*(dist*3)

        %q(1:N/3) = [1:N/3]'*(r_init)
        %q(N/3:2*N/3) = [N/3:2*N/3]'*(r_init+2) + 25
        %q(2*N/3:N) = [2*N/3:N]'*(1) + 25
        
        for c0 = 1:m
            %q_goal(:,c0) = -25 + [0:N-1]'.*r_lattice; %based on node 1 as goal
            goal = -25;
            q_goal(:,c0) = goal + [-leadNode + 1: N - leadNode]'.*r_lattice;
            %q_goal(leadNode,c0) = goal;
            %q_goal(1:leadNode-1,c0) = goal + [0:leadNode-1]'.*r_lattice;
            %q_goal(leadNode+1:N,c0) = goal + [leadNode+1:N]'.*r_lattice;
        end
        leadNode
        N
        q_goal
        q_goal(leadNode)
    end

    %q=[0:(r_lattice-delta):(r_lattice-delta)*(N-1)]';
    
    %q = [5; 50; 90; 120; 150; 180; 210; 230; 250; 260]
    %q = [1:N]'*r_init
    %q = [5; 5+r_safety+((r_init-r_safety)/2); 19; 26; 40; 47; 55; 62; 70; 77]
    %p = [v_max; -v_max; v_max; v_max; v_max; -v_max; v_max; -v_max; -v_max; v_max]
    %p = [v_max; -v_max; v_max; v_max; v_max; -v_max; v_max; -v_max; -v_max; v_max;
    %     v_max; -v_max; v_max; v_max; v_max; v_max; v_max: v_max: v_max; -v_max; -v_max; -v_max]
    %p = [v_max; -v_max];
    %v_max; v_max; v_max; v_max; v_max; v_max; -v_max; v_max; 
        %-v_max; v_max; v_max; v_max; v_max; v_max; v_max; v_max; v_max; -v_max]

    spatial_neighbors = zeros(N, N, N);

    tcyc=Tc;
    time_ctrl=[0:tcyc:tmax-tcyc]';       %vector of control cycle times
    time_traj=[0:tcyc/tdiv:tmax]';  %vector of trajectories
    u_steps=round(((tmax-tcyc)/tcyc))+1;   %number of control update steps
    steps=tmax/(tcyc/tdiv);         %total numbers of trajectory steps

    %goals for gamma control term
    qd=zeros(N, m);         %preallocate
    pd=zeros(N, m); %end with 0 velocity
    
    %initial control value
    %u = (a_max + (-a_max - a_max).*rand(N, m));
    if framework ~= 2
        u = zeros(N, m);
    else
        u = q;
    end
    uGradient = zeros(N, m);
    uConsensus = zeros(N, m);
    uGamma = zeros(N, m);
    uBetaGradient = zeros(N, m);
    uBetaConsensus = zeros(N, m);
    uNew = zeros(N, m);

%     %generate an equally spaced (by r_lattice) rectangular grid goal
%     %TODO: generalize for m-dimensional (use cat and m-th root of N [need m loops?])
%     if m == 1
%         for i=1:N
%             qd(i,:) = ((i-1)*r_lattice) + coord_max;
%         end
%         %scatter(qd(:,1),zeros(N,1),[1:N]');
%     elseif m == 2
%         for i=1:ceil(sqrt(N))
%             for j=1:ceil(sqrt(N))
%                 qd((i-1)*ceil(sqrt(N))+j,:) = [((i-1)*r_lattice) ((j-1)*r_lattice)] + [coord_max coord_max];
%             end
%         end
%         %scatter(qd(:,1),qd(:,2));
%     elseif m == 3
%         for i=1:ceil(N^(1/3))
%             for j=1:ceil(N^(1/3))
%                 for k=1:ceil(N^(1/3))
%                     qd((i-1)*ceil(N^(2/3))+(j-1)*ceil(N^(1/3))+k,:) = [((i-1)*r_lattice) ((j-1)*r_lattice) ((k-1)*r_lattice)] + [coord_max coord_max coord_max];
%                 end
%             end
%         end
%         %scatter(qd(:,1),qd(:,2),qd(:,3));
%     end
%     qd=qd(1:N,:); %shrink to maximum N

    %qd=ones(N,m).*coord_max*(2+N);
    qd=ones(N,m).*coord_max;
    pd=ones(N,m).*(v_max/5);
    %pd=zeros(N,m)
    
    qr=qd;
    pr=pd;
    
    if framework ~= 2
        q_goal = qd;
    end
    
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
    
    de = deviationEnergy(framework, q, q_goal, r_comm, r_lattice, delta)
    de_norm = de / norm(de, 2)
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
    e_history = zeros([round(steps), N, m]);
    for a = 1 : numLyap+1
        v_history = zeros([a round(steps), 1, 1]); %R^(N*m) -> R
        vdot_history = zeros([a round(steps), 1, 1]);
    end
    alp_history = zeros([round(steps), N, m]);
    p_history = zeros([round(steps), N, m]);
    qr_history = zeros([round(steps), N, m]);
    pr_history = zeros([round(steps), N, m]);
    de_history = zeros([round(steps), N, m]);
    de_norm_history = zeros([round(steps), N, m]);
    
    %uPeriod = (1:N)'.*Tc %different update periods for all particles
    uPeriod = ones(N,1)*Tc %same update period for all particles
    %uPeriod(1) = Tc*5; %make one node update its control slowly
    uOffset = zeros(N,1)

    if delay ~= 0
        %uOffset = delay_max + (delay_min - delay_max).*rand(N, 1) %add arbitrary delay to each control update
        tadd = tcyc / tdiv;
        for i=1:N
            to=round(rand(1,1)*tdiv);
            for j=1:to
                uOffset(i) = uOffset(i) + tadd;
            end
        end
        
        for i=1:N
            if uOffset(i) <= 0
                uOffset(i) = Tc
            end
        end
        
        uOffset
    end

    switches = 0; %initialize to 0 (before any goal updates occur)

    %system evolution
    for t=0:tcyc:tmax-tcyc
        t_i=round(t/tcyc)+1;
        
%         if delay == 1 && t < Tc
%             continue;
%         end

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
                %subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),zeros(N,1),'b');
                scatter(q(:,1),zeros(N,1),'g');
                quiver(q(:,1),zeros(N,1),p(:,1),zeros(N,1),'g');
                scatter(qr(:,1),zeros(N,1),'r'); %plot goal
            elseif m == 2
                quiver(q(:,1),q(:,2),p(:,1),p(:,2),'g')
                subplot(subplotRows,subplotCols,subplotCount), scatter(q(:,1),q(:,2),'g');
                quiver(qr(:,1),qr(:,2),pr(:,1),pr(:,2),'r')
                subplot(subplotRows,subplotCols,subplotCount), scatter(qr(:,1),qr(:,2),'r');
                %out=[[qd(:,1) (qd(:,1)+pd(:,1))] [qd(:,2) (qd(:,2)+pd(:,2))]];
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplot(subplotRows,subplotCols,subplotCount), plot([qd(:,1) (qd(:,1)+pd(:,1))],[qd(:,2) (qd(:,2)+pd(:,2))],'g')
                %subplotCount = subplotCount + 1;
            elseif m == 3
                quiver3(q(:,1),q(:,2),q(:,3),p(:,1),p(:,2),p(:,3),'g');
                subplot(subplotRows,subplotCols,subplotCount), scatter3(q(:,1),q(:,2),q(:,3),'g');
                quiver3(qr(:,1),qr(:,2),qr(:,3),pr(:,1),pr(:,2),pr(:,3),'r');
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
            
            et = errorTransform(framework, q, r_lattice, q_goal(leadNode,:), leadNode);
            %et = et(2:N);
            %et = 1000*(errorTransform(framework, q, r_lattice, q_goal(1,:)).^2) + 0.1*((q - q_goal).^2);
            %et = 1000*errorTransform(framework, q, r_lattice, q_goal(1,:)) + 0.1*((q - q_goal));
            %et = 10000*errorTransform(framework, q, r_lattice, q_goal(1,:)) + 0.01*((q - q_goal));
            %e_history(t_j,:,:) = [0; et];
            %e_history(t_j,:,:) = q - q_goal;
            e_history(t_j,:,:) = et;
            
            if framework == 2
                %alp(:,:) = max(alp_min, min(rand(N,m), alp_max)); %cool effect: any alpha between 0 and 1 works
                alp(:,:) = ones(N,m)*0.9;
                %alp(leadNode) = min(alp);
                %A = diag(1 - alp(:,:));
                A = zeros(N);
                %make the diagonal matrix
                for z = 1 : N
                    if z == leadNode
                        A(z, z) = 1 - alp(1)/2;
                        A(z, z) = 0;
                    else
                        if z ~= leadNode && z < N
                            A(z, z) = 1 - alp(z+1)/2 - alp(z)/2;
                        else
                            A(z, z) = 1 - alp(z) - alp(z - 1)/2;
                        end
                        A(z+1, z) = alp(z)/2;
                        A(z, z+1) = alp(z)/2;
                    end
                end
                A = A(1:N, 1:N);

                %expm(A)

                %logm(expm(A))

                %logm(A)

                %rhobar = max(abs(eig(A)));
                rhobar = norm(eig(A), inf);
                if rhobar > 1
                    'A unstable';
                end

                size(A);
                Q = eye(size(A));
                P = dlyap(A', Q);
                if min(eig(P), 0) < 0
                    'P not positive definite'
                end

                %calculate the next time to jump for node 1 (requires fixed alpha k)
                if t_i == 1
                    %jumpError
                    %norm(e_history(1,:,:), inf)
                    %rhobar
                    %A
                    kjump = ceil(log( jumpError / norm(et(2:N), inf))/log(rhobar)); %todo: change to leadnode
                    %return
                elseif t_i > kjump
                    %if norm(et, inf) <= jumpError
                    norm(e_history(t_i,2:N,:), inf);
                    (log( jumpError / norm(e_history(t_i,2:N,:), inf))/log(rhobar));
                    kjump = kjump + ceil(log( jumpError / norm(e_history(t_i,2:N,:), inf))/log(rhobar)); %todoL change to lead node
                    %return
                end
                
                %passino model specific history
                alp_history(t_j,:,:) = alp(:,:);
            elseif framework == 3
                P = 0;
            end
            
            %store all state variables over time
            q_history(t_j,:,:) = q(:,:);

            if framework == 2
                v_history(1, t_j,:,:) = sum(et'*et);
                %v_history(1, t_j,:,:) = sum((q - q_goal).^2);
                %v_history(1, t_j,:,:) = sum(abs(et(1)).*abs(et(2:N).^2));

                v_history(2, t_j,:,:) = et(1:N)'*P*et(1:N); %todo, leadnode
                %v_history(2, t_j,:,:) = q'*P*q;
                %v_history(2, t_j,:,:) = et'*P*et;

                v_history(3, t_j,:,:) = sum(et.^2);
                %v_history(3, t_j,:,:) = (abs(et(1)).^(N + mod(N,2)) + sum(et(2:N).^2))/2; %force to be even power
                %v_history(3, t_j,:,:) = sum(abs(et));
                %v_history(3, t_j,:,:) = ((et(1)).^(2) + sum(et(2:N).^2)/N);

                if t_j <= 1
                    vdot_history(1, t_j,:,:) =  0;
                    vdot_history(2, t_j,:,:) =  0;
                    vdot_history(3, t_j,:,:) =  0;
                else
                    vdot_history(1, t_j,:,:) = v_history(1, t_j, :, :) - v_history(1, t_j - 1, :, :);

                    vdot_history(2, t_j,:,:) = v_history(2, t_j, :, :) - v_history(2, t_j - 1, :, :);
                    %vdot_history(2, t_j,:,:) = q'*(A'*P + P*A)*q;

                    vdot_history(3, t_j,:,:) = v_history(3, t_j, :, :) - v_history(3, t_j - 1, :, :);
                end
            end
            p_history(t_j,:,:) = p(:,:);
            qr_history(t_j,:,:) = qr(:,:);
            pr_history(t_j,:,:) = pr(:,:);
            de_history(t_j,:,:) = de(:,:);
            de_norm_history(t_j,:,:) = de_norm(:,:);

            %compute control (based on state vector, possibly delayed, etc)
            for i=1:N
                if (mod(tt, (uPeriod(i) + uOffset(i))) == 0) % || t_j == 1 %uncomment to let start at 
                                                                           %t=0 instead of t=Tc
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
                    if framework ~= 2
                        u(i,:) = zeros(1,m);
                    else
                        u(i,:) = q(i,:);
                    end
                    uGradient(i,:) = zeros(1,m);
                    uConsensus(i,:) = zeros(1,m);
                    uBetaGradient(i,:) = zeros(1,m);
                    uBetaConsensus(i,:) = zeros(1,m);
                    uGamma(i,:) = zeros(1,m);
                    uNew(i,:) = zeros(1,m);
                    
                    %delayed state messages
                    if delay == 0
                        q_delay = q;
                        p_delay = p;
                        u_delay = u;
                    else
                        if framework == 0
                            q_delay_tmp = q_history(round(t_j - (uPeriod(i) + uOffset(i)*tdiv)/Tc + 1),:,:);
                            p_delay_tmp = p_history(round(t_j - (uPeriod(i) + uOffset(i)*tdiv)/Tc + 1),:,:);
                        elseif framework == 1
                            q_delay_tmp = q_history(t_j - tdiv + 1,:,:);
                            p_delay_tmp = p_history(t_j - tdiv + 1,:,:);
                            u_delay_tmp = u_history(t_i,:,:);
                        elseif framework == 2
                            %round(t_j - ((uPeriod(i) + uOffset(i)*tdiv))/Tc - tdiv + 2)
                            q_delay_tmp = q_history(round(t_j - ((uPeriod(i) + uOffset(i)*tdiv))/Tc - tdiv + 2),:,:);
                            p_delay_tmp = p_history(round(t_j - ((uPeriod(i) + uOffset(i)*tdiv))/Tc - tdiv + 2),:,:);
                            u_delay_tmp = u_history(t_i,:,:);
                        end

                        if m == 1
                            q_delay(1:N,:) = q_delay_tmp';
                            p_delay(1:N,:) = p_delay_tmp';
                            u_delay(1:N,:) = u_delay_tmp';
                        else
                            q_delay(1:N,:) = q_delay_tmp(:,1:N,:);
                            p_delay(1:N,:) = p_delay_tmp(:,1:N,:);
                            u_delay(1:N,:) = u_delay_tmp(:,1:N,:);
                        end
                    end

                    if framework == 0 %olfati-saber control
                        js = neighborsSpatial(i, q(i,:), q, r_comm, r_lattice);
                        %betas = neighborsSpatial(i, q(i,:), q, r_comm_prime, r_lattice_prime);
                        %js = neighborsSpatialLattice(i, q(i,:), q, r_comm, r_lattice, delta);
                        %if size(js,1) > -1
                            for j=1:size(js,1)
                                %TODO: verify equations and all functions
                                %u(i,:) = u(i,:) + phi_a(sig_norm ( q(j,:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * (((q(j,:) - q(i,:)) ) / (sqrt(1 + epsilon * ((norm(q(j,:) - q(i,:), 2))^2))));
                                %u(i,:) = u(i,:) + (a_ij(q(i,:), q(j,:), r_comm_sig, ha,epsilon) * ((p(j,:) - p(i,:))));
                                %old: uGradient(i,:) = uGradient(i,:) + phi_a(sig_norm ( q(js(j),:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * n_ij(q(i,:), q(js(j),:), epsilon);
                                %old: uConsensus(i,:) = uConsensus(i,:) + (a_ij(q(i,:), q(js(j),:), r_comm_sig, ha, epsilon) * ((p(js(j),:) - p(i,:))));
                                uGradient(i,:) = uGradient(i,:) + phi_a(sig_norm ( q_delay(js(j),:) - q_delay(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * n_ij(q_delay(i,:), q_delay(js(j),:), epsilon);
                                uConsensus(i,:) = uConsensus(i,:) + (a_ij(q_delay(i,:), q_delay(js(j),:), r_comm_sig, ha, epsilon) * ((p_delay(js(j),:) - p_delay(i,:))));
                            end

                         %obstacles
                         %   for j=1:size(betas,1)
                         %       uBetaGradient(i,:) = uBetaGradient(i,:) + phi_a(sig_norm ( q_delay(betas(j),:) - q_delay(i,:), epsilon ), r_comm_prime, r_lattice_prime, ha, a, b) * n_ij(q_delay(i,:), q_delay(betas(j),:), epsilon);
                         %       uBetaConsensus(i,:) = uBetaConsensus(i,:) + (a_ij(q_delay(i,:), q_delay(betas(j),:), r_comm_prime, ha, epsilon) * ((p_delay(betas(j),:) - p_delay(i,:))));
                         %   end
            %             else
            %                 %no neighbors yet: have to compute something
            %                 %TODO: fix, not right (always 0 obviously)
            %                 u(i,:) = u(i,:) + phi_a(sig_norm ( q(i,:) - q(i,:), epsilon ), r_comm_sig, r_lattice_sig, ha, a, b) * (((q(i,:) - q(i,:)) ) / (sqrt(1 + epsilon * norm(q(i,:) - q(i,:), 2))^2));
            %                 u(i,:) = u(i,:) + (a_ij(q(i,:), q(i,:), r_comm_sig, ha, epsilon) * ((p(i,:) - p(i,:))));
            %             end

                        %add gamma goal term
                        %u(i,:) = u(i,:) - c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));
                        %old: uGamma(i,:) = -c1*(q(i,:) - qr(i,:)) - c2*(p(i,:) - pr(i,:));

                        %uGamma(i,:) = -c1gamma*sigma_1((q_delay(i,:) - qr(i,:))) - c2gamma*(p_delay(i,:) - pr(i,:));
                        uGamma(i,:) = -c1gamma*((q_delay(i,:) - qr(i,:))) - c2gamma*(p_delay(i,:) - pr(i,:));
                        
                        %sum all forces
                        %u(i,:) = uGradient(i,:) + uConsensus(i,:) + uGamma(i,:) + c1beta*uBetaGradient(i,:) + c2beta*uBetaConsensus(i,:);
                        u(i,:) = uGradient(i,:) + uConsensus(i,:) + uGamma(i,:);
                    elseif framework == 1
                        %1: check all nodes to see if they see anyone on their
                        %   right, then have them go to goal with max
                        %   acceleration
                        %   for implementation, just check THIS node (i) and
                        %   see if this is the case
                        %2: find nearest right hand neightbor and apply 
                        %   u_i = a * (xhat_i - x_i - s) + b * (vhat_i - v_i)
                        %   where s is lattice separation (r_lattice)

                        %uNew(i,:) = ((qd(i,:) - q(i,:)) / norm(qd(i,:) - q(i,:), 2)) * a_max; %not delayed (can read own state at any time)
                        uNew(i,:) = ((qr(i,:) - q_delay(i,:)) / norm(qr(i,:) - q_delay(i,:), 2)) * a_max; %delayed (can only read own state at sends)
                        dist_min_right = inf;
                        dist_min_left = inf;
                        min_j = -1;

                        %overwrite control if necessary
                        for j=1:N
                            if (norm(q_delay(j) - q_delay(i),2) <= r_comm) && ((q_delay(j)) > (q_delay(i))) && (j ~= i) && (norm(q_delay(j) - q_delay(i),2) <= dist_min_right)
                                %uNew(i,:) = aS * (q(j,:) - q(i,:) - r_lattice + delta) + bS * (p(j,:) - p(i,:));
                                uNew(i,:) = ((1/Tc)^2)*(q_delay(j,:) - q_delay(i,:) - r_init) + ((1/Tc))*(p_delay(j,:) - p_delay(i,:));
                                %uNew(i,:) = (q_delay(j,:) - q_delay(i,:) - r_init) + (p_delay(j,:) - p_delay(i,:));
                                %uNew(i,:) = a_max;
                                dist_min_right = norm(q_delay(j) - q_delay(i),2);
                                %if norm(q(j) - q(i),2) <= (r_lattice - delta)
                                if norm(q_delay(j) - q_delay(i),2) <= r_init
                                    uNew(i,:) = -a_max;
                                end
                                min_j = j; %update each time until last it is set properly
    %                         elseif (norm(q_delay(j) - q_delay(i),2) <= r_comm) && ((q_delay(j)) < (q_delay(i))) && (j ~= i) && (norm(q_delay(j) - q_delay(i),2) <= dist_min_left)
    %                             %slow the lead node down to form flocking
    %                             dist_min_left = norm(q_delay(j) - q_delay(i),2);
    %                             %uNew(i,:) = -a_max;
    %                             uNew(i,:) = -((1/Tc)^2)*(q_delay(j,:) - q_delay(i,:) - r_init) - ((1/Tc))*(p_delay(j,:) - p_delay(i,:));
    %                             %uNew(i,:) = -(aS/4) * (q(j,:) - q(i,:) - r_init) - bS * (p(j,:) - p(i,:));
    %                             if norm(q_delay(j) - q_delay(i),2) <= r_init
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
                        
                        if (min_j > 0) && (a_max <= (2*r_safety)/(Tc^2) + (4*v_max^2)/(a_max * Tc^2) + (6*v_max)/Tc + u_delay(min_j))
                            'Insatisfiable condition';
                            (2*r_safety)/(Tc^2) + (4*v_max^2)/(a_max * Tc^2) + (6*v_max)/Tc + u_delay(min_j);
                            (2*r_safety)/(Tc^2);
                            (4*v_max^2)/(a_max * Tc^2);
                            (6*v_max)/Tc;
                            u_delay(min_j);
                        end

                        if constrain ~= 0
                            if (p_delay(i) >= v_max)
                                uNew(i,:) = -a_max;
                            elseif (p_delay(i) <= -v_max)
                                uNew(i,:) = a_max;
                            end
                            
                            u = sign(u).*(min(abs(u),a_max)); %saturate control
                        end
                    elseif framework == 2
                        dist_min_right = inf;
                        dist_min_left = inf;
                        min_j_right = -1;
                        min_j_left = -1;

                        %overwrite control if necessary
                        for j=1:N
                            if (norm(q_delay(j) - q_delay(i),2) <= r_comm) && ((q_delay(j)) > (q_delay(i))) && (j ~= i) && (norm(q_delay(j) - q_delay(i),2) <= dist_min_right)
                                dist_min_right = norm(q_delay(j) - q_delay(i),2);
                                min_j_right = j; %update each time until last it is set properly
                            elseif (norm(q_delay(j) - q_delay(i),2) <= r_comm) && ((q_delay(j)) < (q_delay(i))) && (j ~= i) && (norm(q_delay(j) - q_delay(i),2) <= dist_min_left)
                                dist_min_left = norm(q_delay(j) - q_delay(i),2);
                                min_j_left = j;
                            end
                        end
                        
                        if i == leadNode
                            %max_sep = 0;
                            etl = [et(1:leadNode-1); et(leadNode+1:N)]; %todo: fix, this will error on leadNode=1
                            max_sep = norm(etl, inf); %todo: lead node
                            
%                             for asdf=2:N
%                                 %error = q_delay(asdf) - (q_delay(asdf - 1) + dist);%this is an optimized control, allowing the first node to go ahead
%                                                                                    %and move before StrongFlock, if the flock is just squished, which results
%                                                                                    %in faster convergence;
%                                                                                    %the absolute version forces the flock to expand first, then move, 
%                                                                                    %NOTE: this is necessary to use the sum of squares common Lyapunov function
%                                 error = abs(q_delay(asdf) - (q_delay(asdf - 1) + dist));%resulting in slower convergence
%                                 
%                                 if (error > max_sep)
%                                     max_sep = error;
%                                 end
%                             end
                            %if et(1) <= 2*jumpError && max_sep <= jumpError
                            %    u(i,:) = q_delay(i) - 0.1*(q_delay(i) - q_goal(1,:));
                            %elseif max_sep <= jumpError
                            if max_sep <= jumpError
%                                 max_sep
%                                 jumpError
%                                 kjump
%                                 t_i
%                                 t_j
                                %return
                                
                                flocking_weak_count = flocking_weak_count + 1;
                                
                                %if flocking_weak_count >= 5
                                    %uchange = q_delay(i) - alp(i,:) * (q_delay(i) - q_goal(1,:));
                                    %uchange = q_delay(i) - alp(i,:) * (q_delay(i) - q_goal(1,:));
                                    %uchange = q_delay(i) - (alp(i,:)/N)*(q_delay(i) - q_goal(1,:))

                                    %uchange = max(q_delay(i - 1) + eps, min(q_delay(i) - (alp(i,:))*(q_delay(i) - q_goal(i,:)), q_delay(i + 1) - eps))
                                    uchange = q_delay(i) - (alp(i,:))*(q_delay(i) - q_goal(i,:));

                                    if abs(q_delay(i) - uchange) < jumpError/2
                                        u(i,:) = uchange;
                                    else
                                        u(i,:) = q_delay(i) + sign(uchange) * jumpError/2;
                                    end

                                    flocking_weak_count = 0;
                                %end
                            else
                                u(i,:) = q_delay(i);
                            end
                            
                            %if norm(q_delay(1) - q_delay(2), 2) <= (dist + eps)
                            %    u(i,:) = q_delay(i);
                            %else
                            %    u(i,:) = q_delay(i); %update to same location
                            %end
                        elseif i == N
                            %u(i,:) = q_delay(N) - alp(i,:) * (q_delay(N) - (q_delay(N - 1) + dist));
                            uchange = max(q_delay(N - 1) + eps, q_delay(N) - alp(i,:) * (q_delay(i) - (q_delay(i - 1) + dist)));
                            if abs(q_delay(i) - uchange) < step_max
                                u(i,:) = uchange;
                            else
                                u(i,:) = q_delay(i) + sign(uchange) * step_max;
                            end
                        elseif i == 1
                            uchange = min(q_delay(i + 1) - eps, q_delay(i) - alp(i,:) * (q_delay(i) - (q_delay(i + 1) - dist)));
                            if abs(q_delay(i) - uchange) < step_max
                                u(i,:) = uchange;
                            else
                                u(i,:) = q_delay(i) + sign(uchange) * step_max;
                            end
                        else %2...N-1
                            %u(i,:) = q_delay(i) - alp(i,:) * (q_delay(i) - (1/2)*(q_delay(i - 1) + q_delay(i + 1)));
                            uchange = max(q_delay(i - 1) + eps, min(q_delay(i) - alp(i,:) * (q_delay(i) - (1/2)*(q_delay(i - 1) + q_delay(i + 1))), q_delay(i + 1) - eps));
                            if abs(q_delay(i) - uchange) < step_max
                                u(i,:) = uchange;
                            else
                                u(i,:) = q_delay(i) + sign(uchange) * step_max;
                            end
                        end
                    end
                end
            end
            
            %run system evolution
            
            %saturate the controls (no effect until out of neighbors loop
            %above, so might as well be efficient and do it once for
            %everything)
            if constrain ~= 0
                u = sign(u).*(min(abs(u),a_max));
                
                %saturate the velocity
                %sign(p)
                %(min(abs(p),v_max))
                p = sign(p).*(min(abs(p),v_max));
            end

            %state space form
            %q'=p;     => [q'; p']=[0 1; 0 0]*[q; p] + [0;1]*[u]
            %p'=u;
            
            %solution is:
            %  expm(A*(t-t0))*x0+int(expm(A*(t-t0))*B*u,tau,t0,t)
            %for
            %  syms q p q0 p0 tau t t0 u
            %  x = [q;p]
            %  x0 = [q0; p0]
            %  A = [0 1; 0 0]
            %  B = [0; 1]
            %which expands to
            %  [x]=[q;= [q0 + p0*(t - t0) + u*(t - t0)^2;
            %       p]   p0 + u*(t - t0)]

            if system_type == 0
                q = q + p.*(tcyc/tdiv) + u.*((tcyc/tdiv)^2);
                p = p + u.*(tcyc/tdiv);
            elseif system_type == 1
                %no velocity, next position is directly computed
                q = u;
            end

            if constrain ~= 0
                p = sign(p).*(min(abs(p),v_max));
                u = sign(u).*(min(abs(u),a_max));
            end
            
            de = deviationEnergy(framework, q, q_goal, r_comm, r_lattice, delta);
            de_norm = de / norm(de,2);

            %gamma agent
            frv = fr(qr, pr);
            if constrain ~= 0
                frvs = sign(frv).*(min(abs(frv),a_max));
            else
                frvs=frv;
            end
            
            if system_type == 0
                qr=qr + pr.*(tcyc/tdiv) + frvs.*((tcyc/tdiv)^2);
                pr=pr + frvs.*(tcyc/tdiv);
            elseif system_type == 1
                qr = frvs;
            end
            
            if constrain ~= 0
                pr = sign(pr).*(min(abs(pr),v_max));
            end
            
            %check for safety violations
%             for xyz=1:N
%                 N_safety = neighborsSpatial(xyz, q(xyz), q, r_safety, r_safety);
%                 if (size(N_safety) > 0)
%                     %'Bad safety'
%                     %N_safety
%                 end
%             end
            
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
        if plotControls >= 3
            %plot distance of each node so we can see how their lattice evolves
            figure;
            hold on;
            if m == 1
                for i=1:N
                    if mod(i,2)==0
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'c');
                        plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'r');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'c','LineWidth',r_comm);
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'r','LineWidth',r_safety);
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)+(r_safety/2),'k--');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)-(r_safety/2),'k--');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)+(r_comm/2),'c');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)-(r_comm/2),'c');
                    else
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'m');
                        plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'b');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'m','LineWidth',r_comm);
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1),'b','LineWidth',r_safety);
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)+(r_safety/2),'g--');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)-(r_safety/2),'g--');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)+(r_comm/2),'m');
                        %plot(time_traj(1:size(q_history(:,i,1))),q_history(:,i,1)-(r_comm/2),'m');
                    end
                end
                legend('q1');
                
                figure;
                hold on;
                for i=1:N
                    if mod(i,2)==0
                        plot(time_traj(1:size(q_history(:,i,1))),de_history(:,i,1),'r');
                    else
                        plot(time_traj(1:size(q_history(:,i,1))),de_history(:,i,1),'b');
                    end
                end
                plot(time_traj(1:size(q_history(:,i,1))),de_norm_history(:,1),'k');
                legend('de1');
                
                figure;
                hold on;
                for i=1:N
                    if mod(i,2)==0
                        plot(time_traj(1:size(q_history(:,i,1))),e_history(:,i,1),'r');
                    else
                        plot(time_traj(1:size(q_history(:,i,1))),e_history(:,i,1),'b');
                    end
                end
                
                for a = 1 : numLyap+1
                    v(a)=plot(time_traj(1:size(q_history(:,i,1))),v_history(a,:,1),'c');
                    vd(a)=plot(time_traj(1:size(q_history(:,i,1))),vdot_history(a,:,1),'c.');
                    %v2=plot(time_traj(1:size(q_history(:,i,1))),v2_history(:,1),'k');
                    %v2d=plot(time_traj(1:size(q_history(:,i,1))),v2dot_history(:,1),'k.');
                    %v3=plot(time_traj(1:size(q_history(:,i,1))),v3_history(:,1),'g');
                    %v3d=plot(time_traj(1:size(q_history(:,i,1))),v3dot_history(:,1),'g.');
                    max(vdot_history(a,:,:))
                    vlast(a) = inf;
                    ibad(a) = 1;

                    for asdf = 1 : size(time_traj(1:size(q_history(:,i,1))))
                        if v_history(a,asdf,1) > vlast
                            %'Error: increasing Lyapunov function'
                            %vlast(a)
                            %v_history(a,asdf,1)
                            vbad(a,ibad(a),:) = [asdf - 2, v_history(a,asdf-1,1)]; %asdf - 2 since time is 0 indexed
                            ibad(a) = ibad(a) + 1;
                            vbad(a,ibad(a),:) = [asdf - 1, v_history(a,asdf,1)];
                            ibad(a) = ibad(a) + 1;
                        end
                        vlast(a) = v_history(a,asdf,1);
                    end
                end
                legend([v,vd],'v', 'vd', 'v2', 'v2d', 'v3', 'v3d');
                %legend([v,vd],'v', 'vd', 'v2', 'v2d', 'v3', 'v3d');
                
                %do this second to get all bad points plotted on very top
                for a = 1 : numLyap+1
                    if ibad(a) > 1
                        plot(vbad(a,:,1), vbad(a,:,2), 'xr');
                    end
                end
                
                if framework == 2
                    figure;
                    hold on;
                    for i=1:N
                        if mod(i,2)==0
                            plot(time_traj(1:size(q_history(:,i,1))),alp_history(:,i,1),'r');
                        else
                            plot(time_traj(1:size(q_history(:,i,1))),alp_history(:,i,1),'b');
                        end
                    end
                    legend('alp1');
                end
            elseif m == 2
                for i=1:N
                    if mod(i,2)==0
                        plot(time_traj,(q_history(:,i,1).^2 + q_history(:,i,2).^2).^(1/2),'r');
                    else
                        plot(time_traj,(q_history(:,i,1).^2 + q_history(:,i,2).^2).^(1/2),'b');
                    end
                    legend('q2');
                end
            end
        end

        if plotControls == 2 || plotControls >= 4
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
                    legend('u', 'uNew', 'q', 'p', 'de');
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
                        plot3(time_traj,q_history(:,i,1) - qr(i,1),q_history(:,i,2) - qr(i,2),'r:');
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
    qr
    pd
    pr
    de
    qc=Ave(q)
    pc=Ave(p)
end


%Obstacle notes

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