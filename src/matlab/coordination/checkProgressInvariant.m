function [ out ] = checkProgressInvariant(q_last, q_now, q_goal, delta)
    N = size(q_now, 1);
    count = 0;
    terminate = 0;

    numMoves = sum((q_last' ~= q_now));
    
    'moves'
    (q_last' ~= q_now)
    
    %todo: add check to see if within quantization error of final goal, if so, not moving is fine, and this detects termination
    %abs(q_now - q_goal);
    if (sum(abs(q_now - q_goal) < (delta).*(1:N)') == N) && (numMoves == 0) %all within error range and no moves to make
        terminate = 2;
        'terminated';
    end
    
    %display which node is the only that can move (to see if there is a pattern)
    if (numMoves > 0 && numMoves < N/2)
        q_last' ~= q_now;
    end 

    out = (numMoves > 0) + terminate;
end
