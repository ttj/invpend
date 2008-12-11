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

%checkExtrema constrains an input value (1x1) to a maximum and minimum
%inputs  valIn  - value to be constrained within valMin and valMax (1x1)
%        valMin - minimum value (1x1)
%        valMax - maximum value (1x1)
%outputs valOut - constrained value to output (1x1)
function [valOut] = checkExtrema(valIn, valMin, valMax)
    if (valIn <= valMin)
        valOut = valMin;
    elseif (valIn >= valMax)
        valOut = valMax;
    else
        valOut = valIn;
    end
end
