%To-do: Need to fill the F function, return the binary for left and right elements
function [parent_pop, left_elems, right_elems] = initialize_pop(low_1, high_1, low_2, high_2, init_pop_size,resolution)
    
    ret = check_valid_range(low_1, high_1,low_2, high_2);
    if ret == false
        parent_pop = [];
        left_elems = [];
        right_elems = [];
        return;
    end
    
    % Compute the maximum digits
    left_digits = cal_number_digits(low_1, high_1,resolution);
    right_digits = cal_number_digits(low_2, high_2, resolution);
    
    % Create a random permutation
    left_init_pop_list = randperm(2^left_digits, init_pop_size);
    right_init_pop_list = randperm(2^right_digits, init_pop_size);
    
    %Convert to binary base
    left_pop_list_bin = de2bi(left_init_pop_list);
    right_pop_list_bin = de2bi(right_init_pop_list);
    
    left_idx = 2.^repmat(left_digits - 1:-1:0, init_pop_size, 1);
    right_idx = 2.^repmat(right_digits - 1:-1:0, init_pop_size, 1);
    
    %Convert to decimal base on formula
    left_idx = left_idx .* left_pop_list_bin;
    right_idx = right_idx .* right_pop_list_bin;
    
    %Sum of Matrix rows
    sum_left_idx = sum(left_idx, 2);
    sum_right_idx = sum(right_idx, 2);
    
    %Compute the real value of initialize pop
    real_init_pop_left_lst = low_1 .+ sum_left_idx * ((high_1 - low_1) / (2^left_digits - 1));
    real_init_pop_right_lst = low_2 .+ sum_right_idx * ((high_2 - low_2)/(2^right_digits - 1));
    
    % TO-DO: Compute the final F(x, y) values
    % Insert F function here.
    parent_pop = compute_func(real_init_pop_left_lst, real_init_pop_right_lst, resolution);
    left_elems = left_pop_list_bin;
    right_elems = right_pop_list_bin;
end

%%%%%%%%%%%%%%%%%%%%************************%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% COMPUTATION FUNCTION %%%%%%%%%%%%%%%%%%%%%
function f = compute_func(x, y, resolution)
    %To-do: Temporary change the equation here.
    f = x.^2 - 4*x.*y + 5*y.^2 - 4*y + 3 + 10^resolution;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Checking the valid input arguments
function ret_status = check_valid_range(l1, u1, l2, u2)
    if l1 < 0
        fprintf('[ERROR]The input lower bound of the first number is negative');
        ret_status = false;
        return;
    elseif u1 < 0
        fprintf('[ERROR]The input upper bound of the first number is negative');
        ret_status = false;
        return;
    elseif l2 < 0
        fprintf('[ERROR]The input lower bound of the second number is negative');
        ret_status = false;
        return;
    elseif u2 < 0
        fprintf('[ERROR]The input upper bound of the second number is negative');
        ret_status = false;
        return;
    end
        
    if l1 > u1
        fprintf('[ERROR]The range of the first input was invalid, lower bound bigger than upper bound');
        ret_status = false;
        return;
    elseif l2 > u2
        fprintf('[ERROR]The range of the second input was invalid, lower bound bigger than upper bound');
        ret_status = false;
        return;
    end
    ret_status = true;
end


function max_number_digits = cal_number_digits(lower_bound, upper_bound, resolution)
    if lower_bound >= upper_bound
        error('Wrong input range, lower bound bigger than upper bound');
        max_number_digits = -1;
        return;
    else
        % We using binary number so we need to do the logarithm base 2 here.
        max_number_digits = ceil(log2((upper_bound - lower_bound) * 10^resolution));
    end
end