function [ out ] = checkProgressInvariant(q_last, q_now, q_goal, quant_beta)
    N = size(q_now, 1);
    count = 0;
    terminate = 0;

    %todo: add check to see if within quantization error of final goal, if so, not moving is fine, and this detects termination
    %abs(q_now - q_goal);
    %quant_beta;
    %(quant_beta*(N-1))*(1:N)';
    %sum(abs(q_now - q_goal) < (2*quant_beta*(N-1))*(1:N)');
    if sum(abs(q_now - q_goal) < (2*quant_beta*(N-1)).*(1:N)') == N %all within error range
    %if sum(abs(q_now - q_goal) < (2*quant_beta*(N-1))) == N %all within error range (with varying quant_beta)
        terminate = 1;
        'terminated';
    end
    out = (sum((q_last' ~= q_now)) > 0) || terminate;
end
