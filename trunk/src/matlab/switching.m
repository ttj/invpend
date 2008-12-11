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

%switching 
%
%inputs current_controller - present controller index (1=safety, 2=base,
%                            3=exp) (1x1)
%       P - P matrix from LMI problem definining stability region, where
%           the matrix is indexed in the 3rd dimension by the current
%           controller (4x4x3)
%       x - current, future, and next future state values (
function [controller] = switching(current_controller, P, x)
    %region based on x'*P*x<1: have to check current controller to determine
    %proper P to use    
    if (x'*P(1:4,1:4,current_controller)*x > 1)
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