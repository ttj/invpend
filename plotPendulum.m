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

%plotPendulum
%inputs titleStr - title string
%       timeVec  - time vector to plot everything against (1xN)
%       xMatrix  - system evolution matrix (4xN)
%       lyapVec  - system Lyapunov function vector (1xN)
%       lyapdotVec - system Lyapunov function derivative vector (1xN)
%       dlyapVec - system discrete Lyapunov function vector (1xN)
%       dlyapdotVec - system discrete Lyapunov function difference equation
%                  vector (1xN)
%       uVecSat  - saturated control vector, that is, constrained to the 
%                  control input's physical constraints, such as the 
%                  maximum and minimum voltage output from a voltage source
%                  (1xN)
%       uVec     - unsaturated control vector (1xN)
function plotPendulum(titleStr, timeVec, xMatrix, lyapVec, lyapdotVec, dlyapVec, dlyapdotVec, uVecSat, uVec)
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

    figure;
    hold on;
    xi=1;
    thetai=2;
    xdoti=3;
    thetadoti=4;
    plot3(xMatrix(xi:xi,:),xMatrix(xdoti:xdoti,:),timeVec, 'b');
    plot3(xMatrix(thetai:thetai,:),xMatrix(thetadoti:thetadoti,:),timeVec, 'r');
    
    
    %stabilizable region
    
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
    
    
    
    tl=length(timeVec);
    for i=1:1:tl
        c1(i)=x_min;
        c2(i)=x_max;
        c3(i)=theta_min_rad;
        c4(i)=theta_max_rad;
        c5(i)=xdot_min;
        c6(i)=xdot_max;
    end;

%    plot3(c1,c3,timeVec);
%    plot3(c2,c4,timeVec);
%    plot3(c1,c5,timeVec);
%    plot3(c2,c6,timeVec);

%     patch([c1 c2],[c3 c4],[timeVec timeVec],[0]) %basic idea of what we want
     %also see: surf([3 -3;3 -3; 3 -3; 3 -3; 3 -3],[-3 -3;-3 -3;-3 -3;-3 -3;-3 -3],[0 0;1 1;2 2;3 3;4 4])

%    surf(c1, c2, timeVec);
%    surf([c1 c2],[c3 c4],[timeVec timeVec])

%grid proper
%[sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);
%plot3(sx(:),sy(:),sz(:),'*r');
%axis(volumebounds(x,y,z,u,v,w))
%grid; box; daspect([2 2 1])
%[sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);

%    alpha(0.25);

    title(titleStr);
    zlabel('Time');
    xlabel('X, Theta');
    ylabel('Xdot, Thetadot');
    legend('X, Xdot', 'Theta, ThetaDot');
    
    figure;
    hold on;
    plot(timeVec, lyapVec, 'g');
    plot(timeVec, lyapdotVec, 'k');
    plot(timeVec, dlyapVec, 'b');
    plot(timeVec, dlyapdotVec, 'c');
    plot(timeVec, xMatrix(xi:xi,:), 'b');
    plot(timeVec, xMatrix(xdoti:xdoti,:), 'g');
    plot(timeVec, xMatrix(thetai:thetai,:), 'b');
    plot(timeVec, xMatrix(thetadoti:thetadoti,:), 'y');
    plot(timeVec, uVecSat, 'r');
    plot(timeVec, uVec, 'r.');
    [badMin, badMax] = checkStateExtrema(xMatrix(1:3,:), [-0.2 deg2rad(-15) -1.0], [0.2 deg2rad(15) 1.0])
    title(titleStr);
    xlabel('Time');
    ylabel('Magnitude');
    legend('V','Vdot','VD','VDdot','x','xdot','theta','thetaDot','u_{sat}','u'); 

end
