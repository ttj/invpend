function [ out ] = checkNodeTypeInvariant(nodeType_last, nodeType_now)
    HEAD_LEADER_NO_FOLLOWERS = 1;
    HEAD_LEADER = 2;
    MIDDLE = 3;
    TAIL_LEADER = 4;

    N = size(nodeType_last', 1);
    count_head_last = 0;
    count_head_now = 0;
    count_middle = 0;
    count_tail_last = 0;
    count_tail_now = 0;
    
    for i = 1 : N
        if nodeType_last(i) == HEAD_LEADER_NO_FOLLOWERS || nodeType_last(i) == HEAD_LEADER
            count_head_last = count_head_last + 1;
        elseif nodeType_last(i) == TAIL_LEADER
            count_tail_last = count_tail_last + 1;
        end
        
        %no new head nodes or tail nodes can be created
        if nodeType_now(i) == HEAD_LEADER_NO_FOLLOWERS || nodeType_now(i) == HEAD_LEADER
            count_head_now = count_head_now + 1;
        elseif nodeType_last(i) == TAIL_LEADER
            count_tail_now = count_tail_now + 1;
        end
            
        if nodeType_last(i) == MIDDLE && nodeType_now(i) ~= MIDDLE %middle nodes are invariant
            count_middle = count_middle + 1;
        end
    end
    
    out = (count_head_last >= count_head_now) && (count_middle == 0) && (count_tail_last >= count_tail_now);
end
