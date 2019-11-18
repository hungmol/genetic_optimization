%To-do: Need to fill the F function
function parent_pop = initialize_pop(first_lower_bound, first_upper_bound,
                        second_lower_bound, second_upper_bound, init_pop_size)
    if nargin < 5 && nargin >=4 
        warning('[WARNING] Not determine the pop size');
        fprintf('The default pop size will be: 100\n');
        init_pop_size = 100;
    elseif nargin < 4
        fprintf('[ERROR]Not enough input arguments\n');
    end
    
    ret = check_valid_range(first_lower_bound, first_upper_bound,second_lower_bound, second_upper_bound);
    if ret == false
        parent_pop = [];
        return;
    end
    
    % Create a random permutation
    left_digits = cal_number_digits(first_lower_bound, first_upper_bound)
    right_digits = cal_number_digits(second_lower_bound, second_upper_bound)
    left_init_pop_list = randperm(2^left_digits, init_pop_size);
    right_init_pop_list = randperm(2^right_digits, init_pop_size);
end

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


function max_number_digits = cal_number_digits(lower_bound, upper_bound)
    if lower_bound >= upper_bound
        error('Wrong input range, lower bound bigger than upper bound');
        max_number_digits = -1;
        return;
    else
        % We using binary number so we need to do the logarithm base 2 here.
        max_number_digits = ceil(log2(10^(upper_bound - lower_bound)));
    end
end
