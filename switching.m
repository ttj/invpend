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

function [controller] = switching(current_controller, Ps, Pb, Pe, x)
    %region based on x'*P*x<1: have to check current controller to determine
    %proper P to use
    if (current_controller == 1)
        P=Ps;
    elseif (current_controller == 2)
        P=Pb;
    else
        P=Pe;
    end;
    
    if (x'*P*x > 1)
       %current controller leaving stabilizable region, switch to Pb or Ps
        if (current_controller == 3)
            controller = 2;
        else
            controller = 1;
        end;
    else
        controller = current_controller;
    end;
end