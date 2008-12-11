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
function plotPendulum(titleStr, timeVec, xMatrix, lyapVec, lyapdotVec, dlyapVec, dlyapdotVec, uVecSat, uVec, xbound, xbound2)
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
    xindex=1;
    thetai=2;
    xdoti=3;
    thetadoti=4;
    plot3(xMatrix(xindex:xindex,:),xMatrix(xdoti:xdoti,:),timeVec, 'b');
    plot3(xMatrix(thetai:thetai,:),xMatrix(thetadoti:thetadoti,:),timeVec, 'r');
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
    plot(timeVec, xMatrix(xindex:xindex,:), 'b');
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
