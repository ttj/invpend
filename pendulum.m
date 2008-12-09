% static float a[4][4] = { {37.62, 58.22, 17.87, 11.61},
% 	{58.22, 313.16, 69.36, 56.09},
% 	{17.87, 69.36, 29.81, 14.81},
% 	{11.61, 56.09, 14.81, 12.04}
% 	};

% xa[0] = 1.0 * (x) + (-0.00051281) * (-theta) + 0.017961 * (xdot) +
% 		(-0.0000026781) * (-thetadot);
% 	xa[1] = 0 * (x) + 1.0056 * (-theta) + 0.0046419 * (xdot) +
% 		0.020029 * (-thetadot);
% 	xa[2] = 0 * (x) + (-0.049519) * (-theta) + 0.80322 * (xdot) +
% 		(-0.00043546) * (-thetadot);
% 	xa[3] = 0 * (x) + 0.55967 * (-theta) + 0.44824 * (xdot) +
% 		1.0048 * (-thetadot);
% 
% 	xa[0] += 0.0003618 * local_volts;
% 	xa[1] += (-0.00082708) * local_volts;
% 	xa[2] += 0.034913 * local_volts;
% 	xa[3] += (-0.079879) * local_volts;


%Inverted Pendulum Analysis
clear all;

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
x01=-0.19; x02=-0.2418; x03=0.3; x04=1.41;
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
Klmi=[7.2671   30.1022   11.7542    5.9111];
KS=Klmi;
KS=KE; %%%%%%%%%%%%%%%%%%%%%%%
det(Plmi); %must be greater than 0 (shows P is positive definite)

%Xdot=Abar*X
Abar3=A+B*KS;
Abar3s=As+Bs*KSs;
Abar2=A+B*KB;
Abar1=A+B*KE;

det(Abar3'*Plmi+Plmi*Abar3); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)
det(Abar2'*Plmi+Plmi*Abar2); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)
det(Abar1'*Plmi+Plmi*Abar1); %must be less than 0 (shows positive negative, and shows lyapunov derivative is negative always)

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
%Abarexpm=expm(Abar3*t);

eigAbar3=eig(Abar3)
eigAbar2=eig(Abar2)
eigAbar1=eig(Abar1)

%Lyapunov and Discrete Lyapunov Equations for each System
P3=lyap(Abar3',Q)
PD3=dlyap(Abar3',Q)
P2=lyap(Abar2',Q)
PD2=dlyap(Abar2',Q)
P1=lyap(Abar1',Q)
PD1=dlyap(Abar1',Q)

SPtBad=0;
SPdtBad=0;

%set up for previous state (for switching)
Sxtmp=X0;
Sutmp=0;

Bxtmp=X0;
Butmp=0;

Extmp=X0;
Eutmp=0;

SwRxtmp=X0;
SwRutmp=0;

fSwRxtmp=X0;
fSwRutmp=0;

%Timing is crucial to our simulation here.  We use the given control cycle
%time of Tc=20ms to control the main operation (that is, the control is 
%updated every 20ms), then tdiv is the number of
%subsamples to also model within this main period, so that we have smoother
%trajectories within this larger control cycle.
tcyc=Tc;
tdiv=30; %minimum is 1 (1 division per control cycle, same as only looking at control cycle)
tmax=5;
time_ctrl=[0:tcyc:tmax]';      %vector of control cycle times
time_traj=[0:tcyc/tdiv:tmax]'; %vector of trajectories

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

Sxtmpi=Sxtmp;

cS=norm(P3*B,inf)/(min(eig(Q)));
dS=norm(KS*A,inf)+norm(KS*B*KS,inf);
tauS=1/(cS*dS);

cB=norm(P2*B,inf)/(min(eig(Q)));
dB=norm(KB*A,inf)+norm(KB*B*KB,inf);
tauB=1/(cB*dB);

cE=norm(P1*B,inf)/(min(eig(Q)));
dE=norm(KE*A,inf)+norm(KE*B*KE,inf);
tauE=1/(cE*dE);

Sxtmpp=[0 0 0 0]'; %starting
Sxtmppp=[0 0 0 0]';

St(1:4,1)=X0;
Bt(1:4,1)=X0;
Et(1:4,1)=X0;

for t=0:tcyc:tmax-tcyc
    i=round(t/tcyc)+1;

    Sxtmpppp=Sxtmppp;
    Sxtmppp=Sxtmpp;
    Sxtmpp=Sxtmp;
    if (i>1)
        SxtmppAprx=[Sxtmpp(1) Sxtmpp(2) ((Sxtmpp(1)-Sxtmpppp(1))/(2*Ts)) ((Sxtmpp(2)-Sxtmpppp(2))/(2*Ts))]';
    else
        SxtmppAprx=Sxtmp; %just assume perfect for starting case (otherwise could use other, but it will cause a step input to system)
    end;
    Sutmp=KS*SxtmppAprx;
    if (i>1)
        Sutmp=KS*Sxtmpp;
    %    Sxtmp=expm(A*t)*Sxtmp + expm(A*t)*B*Sutmp;
    %    Sxtmp=(expm(Abar3*(t))*Sxtmpp); %need to use previous for I.C. (don't update for inside divisions!)
    %    Sxtmp=(A^(t/tcyc))*Sxtmpp + B*Sutmp; %need to use previous for I.C. (don't update for inside divisions!)
        Sxtmp=F*Sxtmpp + G*B*SuSattmp;
    else
        Sutmp=0;
        Sxtmp=X0;
    end;
    SuSattmp=checkExtrema(Sutmp, va_min, va_max); %u=kx: maybe some problem here since we're doing u=k*expm(A+Bk)*x0?
    Sxtmpi=Sxtmp;
    
    if (Sxtmpi'*Plmi*Sxtmpi > 1) %out of safety region
        'bad state:'
        Sxtmpi
    end;
    
    Bxtmpp=Bxtmp;
    Butmp=KB*Bxtmpp;
    BuSattmp=checkExtrema(Butmp, va_min, va_max); %u=kx: maybe some problem here since we're doing u=k*expm(A+Bk)*x0?
    Bxtmp=F*Bxtmpp + G*B*BuSattmp; %need to use previous for I.C. (don't update for inside divisions!)
    Bxtmpi=Bxtmp;
    
    Extmpp=Extmp;
    Eutmp=KE*Extmpp;
    EuSattmp=checkExtrema(Eutmp, va_min, va_max); %u=kx: maybe some problem here since we're doing u=k*expm(A+Bk)*x0?
    Extmp=F*Extmpp + G*B*EuSattmp; %need to use previous for I.C. (don't update for inside divisions!)
    Extmpi=Extmp;
    
    SwRxtmpp=SwRxtmp;
    %Switch between Safety -> Baseline -> Experimental -> Safety properly
    if (mod(i-1,3)==1)
%        SwRxtmp=(expm(Abar3*t)*SwRxtmpp);
%        SwRutmp=checkExtrema(KS*SwRxtmp, va_min, va_max);
        SwRuSattmp=KS*SwRxtmpp;
        SwRutmp=checkExtrema(SwRuSattmp, va_min, va_max);
        SwRxtmp=F*SwRxtmpp + G*B*SwRuSattmp;
    elseif (mod(i-1,3)==2)
%        SwRxtmp=(expm(Abar2*t)*SwRxtmpp);
%        SwRutmp=checkExtrema(KB*SwRxtmp, va_min, va_max);
        SwRuSattmp=KB*SwRxtmpp;
        SwRutmp=checkExtrema(SwRuSattmp, va_min, va_max);
        SwRxtmp=F*SwRxtmpp + G*B*SwRuSattmp;
    else
%        SwRxtmp=(expm(Abar1*t)*SwRxtmpp);
%        SwRutmp=checkExtrema(KE*SwRxtmp, va_min, va_max);
        SwRuSattmp=KE*SwRxtmpp;
        SwRutmp=checkExtrema(SwRuSattmp, va_min, va_max);
        SwRxtmp=F*SwRxtmpp + G*B*SwRuSattmp;
    end;
    SwRxtmpi=SwRxtmp;
    
    for j=(i-1)*tdiv+2 : 1 : tdiv+(i-1)*tdiv+1
        tt=j*(tcyc/tdiv);
%        Sxtmpi=A*Sxtmpi + B*Sutmp;
%        Sxtmp=(expm(Abar3*tt)*X0);
%        Sxtmpi=real((F+G*B*KS)^(tt))*Sxtmp;
%        Sxtmpi=(F+G*B*KS)*Sxtmpi;

        Sxtmpi=Ft*Sxtmpi+Gt*B*SuSattmp;

%        Sxtmpi=(expm(Abar3*tt)*Sxtmp); %need to use previous for I.C. (don't update for inside divisions!)
%        Sxtmpi=(A^(tt))*Sxtmpi + B*Sutmp;
%        Sxtmpi=F*Sxtmpi + G*B*SuSattmp;
%        Sxtmpi=Sxtmp;
        St(1:4,j)=Sxtmpi;
        Sut(j)=Sutmp;     %u
        SuSatt(j)=SuSattmp;
        %SPt(j)=Sxtmpi.'*PF*Sxtmpi; %V=x'Px
        %SDPdt(j)=(Sxtmpi.'*F.'*PF*Sxtmpi) + (Sxtmpi.'*PF*F*Sxtmpi); %Vdot=x'A'Px+x'PAx
        SPt(j)=Sxtmpi.'*P3*Sxtmpi; %V=x'Px
        SPdt(j)=(Sxtmpi.'*Abar3.'*P3*Sxtmpi) + (Sxtmpi.'*P3*Abar3*Sxtmpi); %Vdot=x'A'Px+x'PAx
%        if (j==1)
%            SPdt(j)=0;
%        else
%            SPdt(j)=(SPt(j)-SPt(j-1))/(tcyc/tdiv);
%        end;
        SDPt(j)=Sxtmpi.'*PD3*Sxtmpi; %Vd=x'Px
        SDPdt(j)=(Sxtmpi.'*Abar3.'*PD3*Sxtmpi) + (Sxtmpi.'*PD3*Abar3*Sxtmpi); %Vdot=x'A'Px+x'PAx
        
        %Sxtmp=(expm(Abar2*tt)*X0);
%        Bxtmpi=(expm(Abar2*tt)*Bxtmpi); %need to use previous for I.C. (don't update for inside divisions!)
        Bxtmpi=Ft*Bxtmpi+Gt*B*BuSattmp;
        Bt(1:4,j)=Bxtmpi;
        But(j)=Butmp;     %u
        BuSatt(j)=BuSattmp;
        BPt(j)=Bxtmpi.'*P2*Bxtmpi; %V=x'Px
        BPdt(j)=(Bxtmpi.'*Abar2.'*P2*Bxtmpi) + (Bxtmpi.'*P2*Abar2*Bxtmpi); %Vdot=x'A'Px+x'PAx
        BDPt(j)=Bxtmpi.'*PD2*Bxtmpi; %Vd=x'Px
        BDPdt(j)=(Bxtmpi.'*Abar2.'*PD2*Bxtmpi) + (Bxtmpi.'*PD2*Abar2*Bxtmpi); %Vdot=x'A'Px+x'PAx

        %Extmp=(expm(Abar1*tt)*X0);
%        Extmpi=(expm(Abar1*tt)*Extmpi); %need to use previous for I.C. (don't update for inside divisions!)
%        Extmpi=(F+G*B*KE)*Extmpp;
        Extmpi=Ft*Extmpi+Gt*B*EuSattmp;
        Et(1:4,j)=Extmpi;
        Eut(j)=Eutmp;     %u
        EuSatt(j)=EuSattmp;
        EPt(j)=Extmpi.'*P1*Extmpi; %V=x'Px
        EPdt(j)=(Extmpi.'*Abar1.'*P1*Extmpi) + (Extmpi.'*P1*Abar1*Extmpi); %Vdot=x'A'Px+x'PAx
        EDPt(j)=Extmpi.'*PD1*Extmpi; %Vd=x'Px
        EDPdt(j)=(Extmpi.'*Abar1.'*PD1*Extmpi) + (Extmpi.'*PD1*Abar1*Extmpi); %Vdot=x'A'Px+x'PAx

        SwRxtmpi=Ft*SwRxtmpi+Gt*B*SwRuSattmp;
        if (mod(i-1,3)==1)
%            SwRxtmpi=(expm(Abar3*tt)*SwRxtmpi);
%            SwRxtmpi=(F+G*B*KS)*SwRxtmpp;
            SwRPt(j) =SwRxtmpi.'*P3*SwRxtmpi; %V=x'Px
            SwRPdt(j)=(SwRxtmpi.'*Abar3.'*P3*SwRxtmpi) + (SwRxtmpi.'*P3*Abar3*SwRxtmpi); %Vdot=x'A'Px+x'PAx
        elseif (mod(i-1,3)==2)
%            SwRxtmpi=(expm(Abar2*tt)*SwRxtmpi);
%            SwRxtmpi=(F+G*B*KB)*SwRxtmpp;
            SwRPt(j) =SwRxtmpi.'*P2*SwRxtmpi; %V=x'Px
            SwRPdt(j)=(SwRxtmpi.'*Abar2.'*P2*SwRxtmpi) + (SwRxtmpi.'*P2*Abar2*SwRxtmpi); %Vdot=x'A'Px+x'PAx
        else
%            SwRxtmpi=(expm(Abar1*tt)*SwRxtmpi);
%            SwRxtmpi=(F+G*B*KE)*SwRxtmpp;
            SwRPt(j) =SwRxtmpi.'*P1*SwRxtmpi; %V=x'Px
            SwRPdt(j)=(SwRxtmpi.'*Abar1.'*P1*SwRxtmpi) + (SwRxtmpi.'*P1*Abar1*SwRxtmpi); %Vdot=x'A'Px+x'PAx
        end;
        Swt(1:4,j)=SwRxtmpi;
        SwRut(j)=SwRutmp;    %u
        SwRuSatt(j)=SwRuSattmp;
    end;
end;

plotPendulum('Safety Controller System Trajectory', time_traj, St, SPt, SPdt, SDPt, SDPdt, SuSatt, Sut);
%plotPendulum('Baseline Controller System Trajectory', time_traj, Bt, BPt, BPdt, SDPt, SDPdt, BuSatt, But);
%plotPendulum('Experimental Controller System Trajectory', time_traj, Et, EPt, EPdt, SDPt, SDPdt, EuSatt, Eut);
%plotPendulum('Switching Proper (S->B->E) Controller System Trajectory', time_traj, Swt, SwRPt, SwRPdt, SDPt, SDPdt, SwRuSatt, SwRut);

%Abar3sLyap=Abar3s.'*P+P*Abar3s+Q
%Abar3try=A+B*KSs
%Abar3sLyap=Abar3try.'*P+P*Abar3try+Q

%sol=solve(strcat(char(Abar3sLyap(1,1)),'=0'),strcat(char(Abar3sLyap(1,2)),'=0'),strcat(char(Abar3sLyap(1,3)),'=0'),strcat(char(Abar3sLyap(1,4)),'=0'),strcat(char(Abar3sLyap(2,1)),'=0'),strcat(char(Abar3sLyap(2,2)),'=0'),strcat(char(Abar3sLyap(2,3)),'=0'),strcat(char(Abar3sLyap(2,4)),'=0'),strcat(char(Abar3sLyap(3,1)),'=0'),strcat(char(Abar3sLyap(3,2)),'=0'),strcat(char(Abar3sLyap(3,3)),'=0'),strcat(char(Abar3sLyap(3,4)),'=0'),strcat(char(Abar3sLyap(4,1)),'=0'),strcat(char(Abar3sLyap(4,2)),'=0'),strcat(char(Abar3sLyap(4,3)),'=0'),strcat(char(Abar3sLyap(4,4)),'=0'), 'p11', 'p12', 'p13', 'p14', 'p22', 'p23', 'p24', 'p33', 'p34', 'p44')
%sol=solve(Abar3sLyap(1,1),Abar3sLyap(1,2),Abar3sLyap(1,3),Abar3sLyap(1,4),Abar3sLyap(2,1),Abar3sLyap(2,2),Abar3sLyap(2,3),Abar3sLyap(2,4),Abar3sLyap(3,1),Abar3sLyap(3,2),Abar3sLyap(3,3),Abar3sLyap(3,4),Abar3sLyap(4,1),Abar3sLyap(4,2),Abar3sLyap(4,3),Abar3sLyap(4,4))



%stabilized by switching but not others?
% A =
% 
%                    0   1.000000000000000                   0                   0
%                    0  -0.100000000000000  -0.010000000000000   0.010000000000000
%                    0                   0                   0   1.000000000000000
%                    0   0.010000000000000   0.500000000000000  -0.100000000000000
% 
% B
% 
% B =
% 
%                    0
%    0.072874493927126
%                    0
%   -0.109311740890688














%lmi

%clear all;

addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/extras'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/demos'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules/parametric'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules/moment'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules/global'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules/robust'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/modules/sos'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/operators'));
addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/solvers/'));

%need:
%Installmex (for solver)
%startup (for solver)

%uncomment on first run to install and start solver
%old=pwd();
%cd 'C:\Users\tjohnson\Documents\Research\inverted_pendulum\matlab\SDPT3-4.0-beta'
%Installmex;
%startup;
%cd old;

addpath(genpath('C:/Program Files (x86)/MATLAB/R2008a/toolbox/yalmip/solvers/yalmip2sdpa'));
%clear classes;

n=4;
r=1;

A=[0         0    1.0000         0; 0         0         0    1.0000; 0   -2.7500  -10.9500    0.0043; 0   28.5800   24.9200   -0.0440];
B=[0;0;1.9400;-4.4400];

Almi=A;
Blmi=B;

thetaMin=deg2rad(-15);
thetaMax=deg2rad(15);
xMin=-0.2;
xMax=0.2;
xDotMin=-1;
xDotMax=1;
VaMin=-4.95;
VaMax=4.95;

a1=[1/xMax 0 0 0]';
a2=[1/xMin 0 0 0]';
%a3=[0 1/xDotMax 0 0]';
%a4=[0 1/xDotMin 0 0]';
a3=[0 1/thetaMax 0 0]';
a4=[0 1/thetaMin 0 0]';
%a5=[0 0 1/thetaMax 0]';
%a6=[0 0 1/thetaMin 0]';
a5=[0 0 1/xDotMax 0]';
a6=[0 0 1/xDotMin 0]';
b1=[1/VaMax];
b2=[1/VaMin];

Qlmi=sdpvar(n,n);
Ylmi=sdpvar(1,n);
Flmi=set(Qlmi > 0);
Flmi=Flmi+set(Qlmi*Almi' + Almi*Qlmi + Ylmi'*Blmi' + Blmi*Ylmi < 0);
Flmi=Flmi+set(a1'*Qlmi*a1 <= 1);
Flmi=Flmi+set(a2'*Qlmi*a2 <= 1);
Flmi=Flmi+set(a3'*Qlmi*a3 <= 1);
Flmi=Flmi+set(a4'*Qlmi*a4 <= 1);
Flmi=Flmi+set(a5'*Qlmi*a5 <= 1);
Flmi=Flmi+set(a6'*Qlmi*a6 <= 1);
Flmi=Flmi+set([eye(r) b1'*Ylmi; Ylmi'*b1 Qlmi] >= 0);
Flmi=Flmi+set([eye(r) b2'*Ylmi; Ylmi'*b2 Qlmi] >= 0);

%sdpsettings('solver','maxdet')
sdpsettings('solver','sdpt3')
sdpsettings('debug',1);

solvesdp(Flmi,-logdet(Qlmi))
Plmi=inv(double(Qlmi))
Klmi=double(Ylmi)*Plmi

a7=[Klmi/VaMax]';
a8=[Klmi/VaMin]';

%now find area
%min logdet(Q^-1)
%s.t.: QA'+AQ<0
%      Q>0
%      aj'Qaj<=1, k=1...6
%A has to be A+BK for QA'+AQ<0
KE=[10.0, 103.36, 27.72, 23.04];
A2lmi=A+B*KE;
Q2lmi=sdpvar(n,n);
F2lmi=set(Q2lmi > 0);
F2lmi=F2lmi+set(Q2lmi*A2lmi'+A2lmi*Q2lmi<0);
F2lmi=F2lmi+set(a1'*Q2lmi*a1 <= 1);
F2lmi=F2lmi+set(a2'*Q2lmi*a2 <= 1);
F2lmi=F2lmi+set(a3'*Q2lmi*a3 <= 1);
F2lmi=F2lmi+set(a4'*Q2lmi*a4 <= 1);
F2lmi=F2lmi+set(a5'*Q2lmi*a5 <= 1);
F2lmi=F2lmi+set(a6'*Q2lmi*a6 <= 1);
F2lmi=F2lmi+set(a7'*Q2lmi*a7 <= 1);
F2lmi=F2lmi+set(a8'*Q2lmi*a8 <= 1);
solvesdp(F2lmi,-logdet(Q2lmi))
P2lmi=inv(double(Q2lmi))

i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=0.01;
xDotStep=0.01;
for (x=xMin:xStep:xMax)
    for (xdot=xDotMin:xDotStep:xDotMax)
        xi=[x 0 xdot 0]';
        if (norm(Klmi*xi) >= VaMax)
            xub(1:4,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb(1:4,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood(1:4,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound(1:4,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

%figure;
%hold on;
%plot(xgood(1:1,:),xgood(3:3,:),'bx')
%plot(xpb(1:1,:),xpb(3:3,:),'g*');
%plot(xub(1:1,:),xub(3:3,:),'r+');
%plot(xbound(1:1,:),xbound(3:3,:),'ko');
i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=0.01;
xDotStep=0.01;
for (x=thetaMin:xStep:thetaMax)
    for (xdot=-2.75:xDotStep:2.75)
        xi=[0 x 0 xdot]';
        if (norm(Klmi*xi) >= VaMax)
            xub2(1:4,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb2(1:4,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood2(1:4,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound2(1:4,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

%figure;
%hold on;
%plot(xgood(2:2,:),xgood(4:4,:),'bx')
%plot(xpb(2:2,:),xpb(4:4,:),'g*');
%plot(xub(2:2,:),xub(4:4,:),'r+');
xbound=xbound(:,2:end);
xbound2=xbound2(:,2:end);
plot(xbound(1:1,:),xbound(3:3,:),'bo'); %x
plot(xbound2(2:2,:),xbound2(4:4,:),'ro'); %theta

%for plotting xbound and xbound2 simultaneously
%cat(2,xbound,zeros(4,length(xbound2)-length(xbound)))

%xboundres1=resample(xbound(1:1,:),length(xbound2),length(xbound),100);
%xboundres3=resample(xbound(3:3,:),length(xbound2),length(xbound),100);
%plot3(xboundres1,xbound2(2:2,:),xbound2(4:4,:),'go');
%plot3(xboundres3,xbound2(2:2,:),xbound2(4:4,:),'yo');


