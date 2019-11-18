function natural_selection(init_pop_f, left_bin_elems, right_bin_elems, pop_size)
    %Compute the probability 
    prob = init_pop_f / sum(init_pop_f, 2);
    
    % Position of probability
    prob_pos = zeros(pop_size,1);
    accum = 0;
    for i = 1:pop_size
        accum += prob(i, 1);
        prob_pos(i,1) = accum;
    end
    

    %%%%%%%%%%%%%%% NEED TO DO THE UNIT TEST FOR THIS FUNCTION %%%%%%%%%%%%
    % Random an array from 0 -> 1
    r_arr = rand(pop_size,1);
    %Init the good pop from initial pop
    better_pop = zeros(pop_size,1);
    %Init for select left and right binary elements.
    better_left_bin = zeros(size(left_bin_elems));
    better_right_bin = zeros(size(right_bin_elems));
    for it=1:pop_size
        r_elem = repmat(r_arr(it,1),pop_size,1);
        r_elem(prob_pos - r_elem <= 0) = nan;
        %finding the minimum bigger than 0
        [r_val, r_index] = min(r_elem,[],1);
        better_pop(it,1) = prob_pos(r_index,1);
        better_left_bin(it,:) = left_bin_elems(r_index,:);
        better_right_bin(it,:) = right_bin_elems(r_index,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

function GEN_crossing(prob_cross, left_bin_elems, right_bin_elems)
    
end