%Inverted Pendulum
%
%Fall 2008
%
%Taylor Johnson
%johnso99@nospam-illinois.edu
%University of Illinois at Urbana-Champaign
%Coordinated Science Laboratory
%Department of Electrical and Computer Engineering

%checkStateExtrema counts the number of elements in a matrix that are
%outside of defined thresholds.  The matrix is of length N and width M and 
%thus composed of M values that are each compared to maxima and minuma.
%inputs  valIn   - vector to compare (MxN)
%        valMin  - minimum values of each element in valIn (Mx1)
%        valMax  - maximum value of each element in valIn (Mx1)
%outputs valMinCount - number of elements in valIn that were below valMin
%                      threshold
%        valMaxCount - number of elements in ValIn that were above valMax
%                      threshold
function [valMinCount, valMaxCount] = checkStateExtrema(valIn, valMin, valMax)
    for j=1 : size(valIn, 1)
        valMinCount(j) = 0;
        valMaxCount(j) = 0;
    end

    for i=1 : size(valIn, 2)
        for j=1 : size(valIn, 1)
            if (valIn(j, i) <= valMin(j))
                valMinCount(j) = valMinCount(j) + 1;
            elseif (valIn(j, i) >= valMax(j))
                valMaxCount(j) = valMaxCount(j) + 1;
            end
        end
    end
end
