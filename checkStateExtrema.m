
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