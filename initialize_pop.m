%To-do: Need to fill the F function
function parent_pop = initialize_pop(first_lower_bound, first_upper_bound,
                        second_lower_bound, second_upper_bound, init_pop_size, resolution)
    if nargin == 4 
        warning('[WARNING] Not determine the pop size and resolution');
        fprintf('The default pop size will be: 100\n');
        fprintf('The default resolution will be: 3\n');
        init_pop_size = 100;
        resolution = 3;
    elseif nargin == 5
        warning('[WARNING] Not determine the resolution');
        fprintf('The default resolution will be: 3\n');
        resolution = 3;
    elseif nargin < 4
        fprintf('[ERROR]Not enough input arguments\n');
        return;
    elseif nargin > 6
        fprintf('Redundance input arguments\n');
        return;
    end
    
    ret = check_valid_range(first_lower_bound, first_upper_bound,second_lower_bound, second_upper_bound);
    if ret == false
        parent_pop = [];
        return;
    end
    
    % Compute the maximum digits
    left_digits = cal_number_digits(first_lower_bound, first_upper_bound,resolution);
    right_digits = cal_number_digits(second_lower_bound, second_upper_bound, resolution);
    
    % Create a random permutation
    left_init_pop_list = randperm(2^left_digits, init_pop_size);
    right_init_pop_list = randperm(2^right_digits, init_pop_size);
    
    %Convert to binary base
    left_pop_list_bin = de2bi(left_init_pop_list);
    right_pop_list_bin = de2bi(right_init_pop_list)
    
    left_idx = 2.^repmat(left_digits - 1:-1:0, init_pop_size, 1);
    right_idx = 2.^repmat(right_digits - 1:-1:0, init_pop_size, 1);
    
    %Convert to decimal base on formula
    left_idx = left_idx .* left_pop_list_bin;
    right_idx = right_idx .* right_pop_list_bin;
    
    %Sum of Matrix rows
    sum_left_idx = sum(left_idx, 2);
    sum_right_idx = sum(right_idx, 2);
    
    %Compute the real value of initialize pop
    real_init_pop_left_lst = first_lower_bound + sum_left_idx * ((first_upper_bound - first_lower_bound) / (2^left_digits - 1));
    real_init_pop_right_lst = second_lower_bound .+ sum_right_idx * ((second_upper_bound - second_lower_bound)/(2^right_digits - 1));
    
    % TO-DO: Compute the final F(x, y) values
end

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
