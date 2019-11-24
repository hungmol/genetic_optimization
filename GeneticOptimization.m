classdef GeneticOptimization
    properties(Access = 'private')
        low_range_1 = 0;
        low_range_2 = 0;
        high_range_1 = 0;
        high_range_2 = 0;
        pop_size = 100; % default
        resolution = 3;
        max_iter = 100;
        
        % Binary and Decimal number parameters
        init_pop_left_bin_lst = [];
        init_pop_right_bin_lst = [];
        left_digits_num = 0;
        right_digits_num = 0;
        
        %Real number in decimal
        init_pop_left_dec_lst = [];
        init_pop_right_dec_lst = [];
        
        % The value of the input function f(x,y)
        f_value_init_lst = [];
        
        % cross breeding
        cross_prob = 20; %unit percentage: 20% by default
    end
    
    methods (Access = 'public')
        %Constructor
        function obj = GeneticOptimization (low_1, high_1, low_2, high_2, pop_size, res, max_iter, cross_prob)
            if(nargin == 8)
                obj.low_range_1     = low_1;
                obj.low_range_2     = low_2;
                obj.high_range_1    = high_1;
                obj.high_range_2    = high_2;
                obj.pop_size        = pop_size;
                obj.resolution      = res;
                obj.max_iter        = max_iter;
                obj.cross_prob      = cross_prob;
            else
                fprintf('[ERROR] Invalid input arguments\n');
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access = 'private')
        %TO-DO: Change the input function equation here
        %%%%%%%%%%%%%%%%%%%%% COMPUTATION FUNCTION %%%%%%%%%%%%%%%%%%%%%
        function f = compute_func(obj, x, y, resolution)
            % Put resolution here to generate the constant C big enough.
            %To-do: Temporary change the equation here.
            f = x.^2 - 4*x.*y + 5*y.^2 - 4*y + 3 + 10^resolution;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = generate_binary_list(obj)
            % Compute the maximum digits
            obj.left_digits_num = cal_number_digits(obj, obj.low_range_1, obj.high_range_1, obj.resolution);
            obj.right_digits_num = cal_number_digits(obj, obj.low_range_2, obj.high_range_2, obj.resolution);
              
            % Create a random permutation
            tmp_init_pop_left_dec_lst = randperm(2^obj.left_digits_num, obj.pop_size);
            tmp_init_pop_right_dec_lst = randperm(2^obj.right_digits_num, obj.pop_size);
                
            %Convert to binary base
            obj.init_pop_left_bin_lst = de2bi(tmp_init_pop_left_dec_lst);
            obj.init_pop_right_bin_lst = de2bi(tmp_init_pop_right_dec_lst);            
        end
        
        function obj = convert_bin_2_dec(obj)
            if obj.left_digits_num - 1 < 0
                fprintf('[ERROR] The number digits of the left side is invalid\n');
                obj = [];
                return;
            elseif obj.right_digits_num - 1 < 0    
                fprintf('[ERROR] The number digits of the right side is invalid\n');
                obj = [];
                return;
            end
            
            left_idx = 2.^repmat(obj.left_digits_num - 1:-1:0, obj.pop_size, 1);
            right_idx = 2.^repmat(obj.right_digits_num - 1:-1:0, obj.pop_size, 1);
            
            %Convert to decimal base on formula
            left_idx = left_idx .* obj.init_pop_left_bin_lst;
            right_idx = right_idx .* obj.init_pop_right_bin_lst;
          
            %Sum of Matrix rows
            sum_left_idx = sum(left_idx, 2);
            sum_right_idx = sum(right_idx, 2);
            
            %Compute the real value of initialize pop
            low_1 = obj.low_range_1;
            low_2 = obj.low_range_2;
            high_1 = obj.high_range_1;
            high_2 = obj.high_range_2;
            left_digits = obj.left_digits_num;
            right_digits = obj.right_digits_num;
            
            obj.init_pop_left_dec_lst = low_1 .+ sum_left_idx * ((high_1 - low_1) / (2^left_digits - 1));
            obj.init_pop_right_dec_lst = low_2 .+ sum_right_idx * ((high_2 - low_2)/(2^right_digits - 1)); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CROSS BREEDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = natural_selection(obj)
            %Compute the probability 
            prob = obj.f_value / sum(obj.f_value, 2); 
            % Position of probability
            prob_pos = zeros(obj.pop_size, 1);
            accum = 0;
            for i = 1:obj.pop_size
                accum += prob(i, 1);
                prob_pos(i,1) = accum;
            end
           
            %%%%%%%%%%%%%%% NEED TO DO THE UNIT TEST FOR THIS FUNCTION %%%%%%%%%%%%
            % Random an array from 0 -> 1
            r_arr = rand(obj.pop_size,1);
            
            %Init the good pop from initial pop
            better_pop = zeros(obj.pop_size,1);
            
            %Init for select left and right binary elements.
            better_pop_left_bin = zeros(size(obj.init_pop_left_bin_lst));
            better_pop_right_bin = zeros(size(obj.init_pop_right_bin_lst));
            for it=1:obj.pop_size
                r_elem = repmat(r_arr(it,1), obj.pop_size, 1);
                r_elem(prob_pos - r_elem < 0) = nan;
                
                %finding the minimum bigger than 0
                [r_val, r_index] = min(r_elem,[],1);
                better_pop(it,1) = prob_pos(r_index,1);
                better_left_bin(it,:) = obj.init_pop_left_bin_lst(r_index,:);
                better_right_bin(it,:) = obj.init_pop_right_bin_lst(r_index,:);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        %%%%%%%%% ULTILITY FUNCTION %%%%%%%%%%%
        %% Checking the valid input arguments
        function ret_status = check_valid_range(obj)
            if obj.low_range_1 < 0
                fprintf('[ERROR]The input lower bound of the first number is negative');
                ret_status = false;
                return;
            elseif obj.high_range_1 < 0
                fprintf('[ERROR]The input upper bound of the first number is negative');
                ret_status = false;
                return;
            elseif obj.low_range_2 < 0
                fprintf('[ERROR]The input lower bound of the second number is negative');
                ret_status = false;
                return;
            elseif obj.high_range_2 < 0
                fprintf('[ERROR]The input upper bound of the second number is negative');
                ret_status = false;
                return;
            end
              
            if obj.low_range_1 >= obj.high_range_1
                fprintf('[ERROR]The range of the first input was invalid, lower bound bigger than upper bound');
                ret_status = false;
                return;
            elseif obj.low_range_2 >= obj.high_range_2
                fprintf('[ERROR]The range of the second input was invalid, lower bound bigger than upper bound');
                ret_status = false;
                return;
            end
            ret_status = true;
         end
         
        function max_number_digits = cal_number_digits(obj, lower_bound, upper_bound, resolution)
            if lower_bound >= upper_bound
                error('Wrong input range, lower bound bigger than upper bound\n');
                max_number_digits = -1;
                return;
            else
                % We using binary number so we need to do the logarithm base 2 here.
                max_number_digits = ceil(log2((upper_bound - lower_bound) * 10^resolution));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUBLIC METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = 'public')
        function run(obj)
            ret = check_valid_range(obj)
            if isequal(ret, false)
                return;
            end
            obj = generate_binary_list(obj);
            obj = convert_bin_2_dec(obj);
            if isempty(obj)
                fprintf('[ERROR] Convert binary to decimal number failed\n');
                return;
            end
            % TO-DO: Compute the final F(x, y) values
            % Insert F function here.
            obj.f_value_init_lst = compute_func(obj, obj.init_pop_left_dec_lst, obj.init_pop_right_dec_lst, obj.resolution);
            disp(obj.f_value_init_lst);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING FIELD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function print_test(obj)
            printf('first low range = %d\n', obj.low_range_1);
            printf('first low range = %d\n', obj.high_range_1);
            printf('second low range = %d\n', obj.low_range_2);
            printf('second high range = %d\n', obj.high_range_2);
            printf('pop size = %d\n', obj.pop_size);
            printf('resolution = %d\n', obj.resolution);
            printf('max iteration = %d\n', obj.max_iter);
        end
        
        function out = test_func(obj, a, b)
            fprintf('Low range 1: %d\n', a);
            fprintf('High range 1: %d\n', b);
            out = a + b;
            fprintf('Output object: %d\n', out);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end 
    
%    methods (Access = 'private')      
%    end
end