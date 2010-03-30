%Inverted Pendulum
%
%Fall 2008
%
%Taylor Johnson
%johnso99@nospam-illinois.edu
%University of Illinois at Urbana-Champaign
%Coordinated Science Laboratory
%Department of Electrical and Computer Engineering
%
%Linear Matrix Inequality
%
%This code utilizes the following third-party modules:
%1. YALMIP Matlab library
%   http://control.ee.ethz.ch/~joloef/wiki/pmwiki.php
%2. SDPT3 Semidefinite Programming (SDP) Solver
%   http://www.math.nus.edu.sg/~mattohkc/sdpt3.html

clear all;

PATH_YALMIP = 'C:/Program Files/MATLAB/R2009b/toolbox/yalmip'; %string path to installation of yalmip
PATH_SDPT3 = 'C:/Program Files/MATLAB/R2009b/toolbox/SDPT3-4.0'; %string path to sdpt3 solver

addpath(genpath(PATH_YALMIP));
addpath(genpath(strcat(PATH_YALMIP, '/extras')));
addpath(genpath(strcat(PATH_YALMIP, '/demos')));
addpath(genpath(strcat(PATH_YALMIP, '/modules')));
addpath(genpath(strcat(PATH_YALMIP, '/modules/parametric')));
addpath(genpath(strcat(PATH_YALMIP, '/modules/moment')));
addpath(genpath(strcat(PATH_YALMIP, '/modules/global')));
addpath(genpath(strcat(PATH_YALMIP, '/modules/robust')));
addpath(genpath(strcat(PATH_YALMIP, '/modules/sos')));
addpath(genpath(strcat(PATH_YALMIP, '/operators')));
addpath(genpath(strcat(PATH_YALMIP, '/solvers/')));

%need:
%Installmex (for solver)
%startup (for solver)

%uncomment on first run to install and start solver
%old=pwd();
%cd PATH_SDPT3
%Installmex;
%startup;
%cd old;

addpath(genpath(strcat(PATH_YALMIP, '/solvers/yalmip2sdpa')));
clear classes;

% %physical parameters
% g=9.8; %gravity: -9.8 m/s
% m=1;    %mass
% Ki=0.9; %torque constant
% Kb=0.9; %back-emf constant
% Kg=50;   %gear ratio
% rad=2;    %
% Ra=0.95;    %armature resistance
% Bm=0.95;    %
% Btheta=0.95;
% Bbar=((Kg*Bm)/(rad^2))+(((Kg^2)*Ki*Kb)/((rad^2)*Ra));
% Bl=((Kg*Ki)/(rad*Ra));
% M=5;    %mass
% Jm=5;   %moment of inertia
% l=0.5;  %length
% Mbar=m+M+(Kg*Jm)/(rad^2); %effective mass
% Dl=(4*Mbar)-3*m;
% 
% a22=(4*Bbar)/(Dl);
% a23=(3*m*g)/(Dl);
% a24=(6*Btheta)/(l*Dl);
% a42=(6*Bbar)/(l*Dl);
% a43=(6*Mbar*g)/(l*Dl);
% a44=(12*Mbar*Btheta)/(m*(l^2)*Dl);
% 
% b2 = 4*Bl/Dl;
% b4 = 6*Bl/(l*Dl);

%A=[0    1.0000 0 0; 0 -10.95 -2.75 0.0043; 0 0 0 1; 0 24.9200 28.5800 -0.0440];
%%A=[0 1 0 0; 0 -a22 -a23 a24; 0 0 0 1; 0 a42 a43 -a44];
%B=[0; 1.94; 0; -4.44];
%%B=[0; b2; 0; -b4]

n=8;
r=2;

a22 = 10.95;
a23 = 2.75;
a24 = 0.0043;
a42 = 24.92;
a43 = 28.58;
a44 = 0.044;

b2 = 1.94;
b4 = 4.44;

A=[0 1.0000 0 0 0 -1 0 0; 0 -a22 -a23 a24 0 0 0 0; 0 0 0 1 0 0 0 0; 0 a42 a43 -a44 0 0 0 0; 0 1 0 0 0 1 0 0; 0 0 0 0 0 -a22 -a23 a24; 0 0 0 0 0 0 0 1; 0 0 0 0 0 a42 a43 -a44]
B=[0 0; b2 0; 0 0; -b4 0; 0 0; 0 b2; 0 0; 0 -b4]

Almi=A;
Blmi=B;

thetaMin=deg2rad(-15);
thetaMax=deg2rad(15);
xMin=-0.2;
xMax=0.2;
xMinB = -0.2;
xMaxB = 0.2;
xDotMin=-1;
xDotMax=1;
VaMin=-4.95;
VaMax=4.95;

VaMaxV = [VaMax; VaMax];
VaMinV = [VaMin; VaMin];

% a1=[1/xMax 0 0 0 1/xMaxB 0 0 0]';
% a2=[1/xMin 0 0 0 1/xMinB 0 0 0]';
% a3=[0 1/xDotMax 0 0 0 1/xDotMax 0 0]';
% a4=[0 1/xDotMin 0 0 0 1/xDotMin 0 0]';
% a5=[0 0 1/thetaMax 0 0 0 1/thetaMax 0]';
% a6=[0 0 1/thetaMin 0 0 0 1/thetaMin 0]';

a1=[1/xMax 0 0 0 0 0 0 0]';
a2=[1/xMin 0 0 0 0 0 0 0]';
a3=[0 1/xDotMax 0 0 0 0 0 0]';
a4=[0 1/xDotMin 0 0 0 0 0 0]';
a5=[0 0 1/thetaMax 0 0 0 0 0]';
a6=[0 0 1/thetaMin 0 0 0 0 0]';
a11=[0 0 0 0 1/xMaxB 0 0 0]';
a12=[0 0 0 0 1/xMinB 0 0 0]';
a13=[0 0 0 0 0 1/xDotMax 0 0]';
a14=[0 0 0 0 0 1/xDotMin 0 0]';
a15=[0 0 0 0 0 0 1/thetaMax 0]';
a16=[0 0 0 0 0 0 1/thetaMin 0]';

%a3=[0 1/thetaMax 0 0]';
%a4=[0 1/thetaMin 0 0]';
%a5=[0 0 1/xDotMax 0]';
%a6=[0 0 1/xDotMin 0]';
%tauS=0.002724511839376;
%tauB=9.946540997091502e-004;
%tauE=6.826194222571933e-004;
%a7=[0 0 0 deg2rad(5)/tauS]';
%a8=[0 0 0 deg2rad(-5)/tauS]';

bc1=[1/VaMax 0]'; %these should be 2x1 vectors, but this doesn't work?
bc2=[1/VaMin 0]';
bc3=[0 1/VaMax]';
bc4=[0 1/VaMin]';

Qlmi=sdpvar(n,n);
Ylmi=sdpvar(r,n);
Flmi=set(Qlmi > 0);
Flmi=Flmi+set(Qlmi*Almi' + Almi*Qlmi + Ylmi'*Blmi' + Blmi*Ylmi < 0);
Flmi=Flmi+set(a1'*Qlmi*a1 <= 1);
Flmi=Flmi+set(a2'*Qlmi*a2 <= 1);
Flmi=Flmi+set(a3'*Qlmi*a3 <= 1);
Flmi=Flmi+set(a4'*Qlmi*a4 <= 1);
Flmi=Flmi+set(a5'*Qlmi*a5 <= 1);
Flmi=Flmi+set(a6'*Qlmi*a6 <= 1);
%Flmi=Flmi+set(a7'*Qlmi*a7 <= 1);
%Flmi=Flmi+set(a8'*Qlmi*a8 <= 1);
Flmi=Flmi+set(a11'*Qlmi*a11 <= 1);
Flmi=Flmi+set(a12'*Qlmi*a12 <= 1);
Flmi=Flmi+set(a13'*Qlmi*a13 <= 1);
Flmi=Flmi+set(a14'*Qlmi*a14 <= 1);
Flmi=Flmi+set(a15'*Qlmi*a15 <= 1);
Flmi=Flmi+set(a16'*Qlmi*a16 <= 1);
Flmi=Flmi+set([eye(1) bc1'*Ylmi; Ylmi'*bc1 Qlmi] >= 0);
Flmi=Flmi+set([eye(1) bc2'*Ylmi; Ylmi'*bc2 Qlmi] >= 0);
Flmi=Flmi+set([eye(1) bc3'*Ylmi; Ylmi'*bc3 Qlmi] >= 0);
Flmi=Flmi+set([eye(1) bc4'*Ylmi; Ylmi'*bc4 Qlmi] >= 0);

%sdpsettings('solver','maxdet')
sdpsettings('solver','sdpt3')
sdpsettings('debug',1);

solvesdp(Flmi,-logdet(Qlmi))
Plmi=inv(double(Qlmi))
Klmi=double(Ylmi)*Plmi


i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=0.005;
xDotStep=0.005;
for x=xMin:xStep:xMax
    for xdot=xDotMin:xDotStep:xDotMax
        xi=[x xdot 0 0 0 0 0 0]';
        if (sum(Klmi*xi >= VaMaxV) > 0 || sum(Klmi*xi <= VaMinV) > 0)
            xub(1:n,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb(1:n,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood(1:n,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound(1:n,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

figure;
hold on;
plot(xgood(1:1,:),xgood(2:2,:),'bx')
plot(xpb(1:1,:),xpb(2:2,:),'g*');
plot(xub(1:1,:),xub(2:2,:),'r+');
plot(xbound(1:1,:),xbound(2:2,:),'ko');
clear xub xpb xgood xbound;



i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=0.005;
xDotStep=0.005;
for x=xMinB:xStep:xMaxB
    for xdot=xDotMin:xDotStep:xDotMax
        xi=[0 0 0 0 x xdot 0 0]';
        if (sum(Klmi*xi >= VaMaxV) > 0 || sum(Klmi*xi <= VaMinV) > 0)
            xub(1:n,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb(1:n,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood(1:n,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound(1:n,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

figure;
hold on;
plot(xgood(5:5,:),xgood(6:6,:),'bx')
plot(xpb(5:5,:),xpb(6:6,:),'g*');
plot(xub(5:5,:),xub(6:6,:),'r+');
plot(xbound(5:5,:),xbound(6:6,:),'ko');
clear xub xpb xgood xbound;


i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=deg2rad(0.5);
xDotStep=deg2rad(0.5);
thetaDotMin=deg2rad(-130);
thetaDotMax=deg2rad(130);
for x=thetaMin:xStep:thetaMax
    for xdot=thetaDotMin:xDotStep:thetaDotMax
        xi=[0 0 x xdot 0 0 0 0]';
        if (sum(Klmi*xi >= VaMaxV) > 0 || sum(Klmi*xi <= VaMinV) > 0)
            xub(1:n,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb(1:n,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood(1:n,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound(1:n,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

figure;
hold on;
plot(rad2deg(xgood(3:3,:)),rad2deg(xgood(4:4,:)),'bx')
plot(rad2deg(xpb(3:3,:)),rad2deg(xpb(4:4,:)),'g*');
plot(rad2deg(xub(3:3,:)),rad2deg(xub(4:4,:)),'r+');
plot(rad2deg(xbound(3:3,:)),rad2deg(xbound(4:4,:)),'ko');
clear xub xpb xgood xbound;





i=1;
iu=1;
ip=1;
ig=1;
ib=1;
bad=0;
badPrev=0;
xStep=deg2rad(0.5);
xDotStep=deg2rad(0.5);
thetaDotMin=deg2rad(-130);
thetaDotMax=deg2rad(130);
for x=thetaMin:xStep:thetaMax
    for xdot=thetaDotMin:xDotStep:thetaDotMax
        xi=[0 0 0 0 0 0 x xdot]';
        if (sum(Klmi*xi >= VaMaxV) > 0 || sum(Klmi*xi <= VaMinV) > 0)
            xub(1:n,iu)=xi;
            iu=iu+1;
            bad=1;
        end;
        if (xi'*Plmi*xi > 1)
            xpb(1:n,ip)=xi;
            ip=ip+1;
            bad=1;
        end;
        if (bad==0)
            xgood(1:n,ig)=xi;
            ig=ig+1;
        end;
        if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
            xbound(1:n,ib)=xi;
            ib=ib+1;
        end;
        i=i+1;
        badPrev=bad;
        bad=0;
    end;
end;

figure;
hold on;
plot(rad2deg(xgood(7:7,:)),rad2deg(xgood(8:8,:)),'bx')
plot(rad2deg(xpb(7:7,:)),rad2deg(xpb(8:8,:)),'g*');
plot(rad2deg(xub(7:7,:)),rad2deg(xub(8:8,:)),'r+');
plot(rad2deg(xbound(7:7,:)),rad2deg(xbound(8:8,:)),'ko');
clear xub xpb xgood xbound;




































% clear a7 a8;
% a7=[0 0 0 deg2rad(45)/tauE]';
% a8=[0 0 0 deg2rad(-45)/tauE]';
% a9=[Klmi/VaMax]';
% a10=[Klmi/VaMin]';
% 
% %now find area
% %min logdet(Q^-1)
% %s.t.: QA'+AQ<0
% %      Q>0
% %      aj'Qaj<=1, k=1...6
% %A has to be A+BK for QA'+AQ<0
% KE=[10.0, 103.36, 27.72, 23.04];
% A2lmi=A+B*KE;
% Q2lmi=sdpvar(n,n);
% F2lmi=set(Q2lmi > 0);
% F2lmi=F2lmi+set(Q2lmi*A2lmi'+A2lmi*Q2lmi<0);
% F2lmi=F2lmi+set(a1'*Q2lmi*a1 <= 1);
% F2lmi=F2lmi+set(a2'*Q2lmi*a2 <= 1);
% F2lmi=F2lmi+set(a3'*Q2lmi*a3 <= 1);
% F2lmi=F2lmi+set(a4'*Q2lmi*a4 <= 1);
% F2lmi=F2lmi+set(a5'*Q2lmi*a5 <= 1);
% F2lmi=F2lmi+set(a6'*Q2lmi*a6 <= 1);
% %F2lmi=F2lmi+set(a7'*Q2lmi*a7 <= 1);
% %F2lmi=F2lmi+set(a8'*Q2lmi*a8 <= 1);
% F2lmi=F2lmi+set(a9'*Q2lmi*a9 <= 1);
% F2lmi=F2lmi+set(a10'*Q2lmi*a10 <= 1);
% solvesdp(F2lmi,-logdet(Q2lmi))
% P2lmi=inv(double(Q2lmi))
% 

% 
% KE=[10.0, 103.36, 27.72, 23.04];
% Abar1=A+B*KE;
% Q=eye(4);
% P1=lyap(Abar1',Q);
% 
% for (x=xMin:xStep:xMax)
%     for (xdot=xDotMin:xDotStep:xDotMax)
%         xi=[x 0 xdot 0]';
%         if (norm(KE*xi) >= VaMax)
%             xub(1:4,iu)=xi;
%             iu=iu+1;
%             bad=1;
%         end;
%         if (xi'*P2lmi*xi > 1)
%             xpb(1:4,ip)=xi;
%             ip=ip+1;
%             bad=1;
%         end;
%         if (bad==0)
%             xgood(1:4,ig)=xi;
%             ig=ig+1;
%         end;
%         if ((bad==0&&badPrev==1) || (bad==1&&badPrev==0))
%             xbound(1:4,ib)=xi;
%             ib=ib+1;
%         end;
%         i=i+1;
%         badPrev=bad;
%         bad=0;
%     end;
% end;
% plot(xbound(1:1,:),xbound(3:3,:),'cs');
% 
