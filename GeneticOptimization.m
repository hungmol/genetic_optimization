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
        
        % New generation list
        new_gen_pop_bin_left_lst = [];
        new_gen_pop_bin_right_lst = [];
        new_gen_val = [];
        
        %Real number in decimal
        init_pop_left_dec_lst = [];
        init_pop_right_dec_lst = [];
        
        % The value of the input function f(x,y)
        f_value_init_lst = [];
        new_f_value_lst = [];
        
        % cross breeding
        cross_prob = 20; %unit percentage: 20% by default
    end
    
    methods (Access = 'public')
        %Constructor
        
        function obj = GeneticOptimization (low_1, high_1, low_2, high_2, pop_size, res, max_iter, cross_prob)
            switch nargin
                case 4
                    obj.low_range_1     = low_1;
                    obj.low_range_2     = low_2;
                    obj.high_range_1    = high_1;
                    obj.high_range_2    = high_2;
                    obj.pop_size        = 100;
                    obj.resolution      = 3;
                    obj.max_iter        = 100;
                    obj.cross_prob      = 20; %percents
                case 5
                    obj.low_range_1     = low_1;
                    obj.low_range_2     = low_2;
                    obj.high_range_1    = high_1;
                    obj.high_range_2    = high_2;
                    obj.pop_size        = pop_size;
                    obj.resolution      = 3;
                    obj.max_iter        = 100;
                    obj.cross_prob      = 20;
                case 6
                    obj.low_range_1     = low_1;
                    obj.low_range_2     = low_2;
                    obj.high_range_1    = high_1;
                    obj.high_range_2    = high_2;
                    obj.pop_size        = pop_size;
                    obj.resolution      = res;
                    obj.max_iter        = 100;
                    obj.cross_prob      = 20;
                case 7
                    obj.low_range_1     = low_1;
                    obj.low_range_2     = low_2;
                    obj.high_range_1    = high_1;
                    obj.high_range_2    = high_2;
                    obj.pop_size        = pop_size;
                    obj.resolution      = res;
                    obj.max_iter        = max_iter;
                    obj.cross_prob      = 20;
                case 8
                    obj.low_range_1     = low_1;
                    obj.low_range_2     = low_2;
                    obj.high_range_1    = high_1;
                    obj.high_range_2    = high_2;
                    obj.pop_size        = pop_size;
                    obj.resolution      = res;
                    obj.max_iter        = max_iter;
                    obj.cross_prob      = cross_prob;
                otherwise
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
            f = x.^2 - 4*x.*y + 5*y.^2 - 4*y + 3 + 10^resolution; % <== this is Constant number
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = generate_binary_list(obj)
            % Compute the maximum digits
            obj.left_digits_num = cal_number_digits(obj, obj.low_range_1, obj.high_range_1, obj.resolution);
            obj.right_digits_num = cal_number_digits(obj, obj.low_range_2, obj.high_range_2, obj.resolution);
              
            % Create a random permutation, size 1 x pop_size
            tmp_init_pop_left_dec_lst = randperm(2^obj.left_digits_num, obj.pop_size);
            tmp_init_pop_right_dec_lst = randperm(2^obj.right_digits_num, obj.pop_size);
                
            %Convert to binary base is matrix have size: pop_size x bit_numbers
            obj.init_pop_left_bin_lst = de2bi(tmp_init_pop_left_dec_lst);
            obj.init_pop_right_bin_lst = de2bi(tmp_init_pop_right_dec_lst);   
        end
        
        function [pop_left_dec, pop_right_dec] = convert_bin_2_dec(obj, pop_left_bin_lst, pop_right_bin_lst)
            if obj.left_digits_num - 1 < 0
                fprintf('[ERROR] The number digits of the left side is invalid\n');
                pop_left_dec = [];
                pop_right_dec = [];
                return ;
            elseif obj.right_digits_num - 1 < 0    
                fprintf('[ERROR] The number digits of the right side is invalid\n');
                pop_left_dec = [];
                pop_right_dec = [];
                return ;
            end
            
            pop_size = size(pop_left_bin_lst)(1,1);
            
            left_idx = 2.^repmat(obj.left_digits_num - 1:-1:0, pop_size, 1);
            right_idx = 2.^repmat(obj.right_digits_num - 1:-1:0, pop_size, 1);
            
            %Convert to decimal base on formula
            left_idx = left_idx .* pop_left_bin_lst;
            right_idx = right_idx .* pop_right_bin_lst;
          
            %Sum of Matrix rows
            sum_left_idx = sum(left_idx, 2);
            sum_right_idx = sum(right_idx, 2);
            
            %Compute the real value of initialize pop
            low_1 = obj.low_range_1;
            low_2 = obj.low_range_2;
            high_1 = obj.high_range_1;
            high_2 = obj.high_range_2;
            
            pop_left_dec = low_1 .+ sum_left_idx * ((high_1 - low_1) / (2^obj.left_digits_num - 1));
            pop_right_dec = low_2 .+ sum_right_idx * ((high_2 - low_2)/(2^obj.right_digits_num - 1)); 
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
        
        % Input is binary number
        function [new_male, new_female] = crossover_a_couple(obj, male, female, crossover_point)
            [num_coup_chromo, num_of_bits] = size(male);
            
            %Male:   1 0 1 0 1 0 | 0 0 0 1 1 1 1 0
            %        ```````````   '''''''''''''''
            %Female: 0 1 0 0 1 0 | 0 1 1 0 0 1 0 1
            %        '''''''''''   ```````````````
            
            % Compute for new male
            male_mask = ones(num_coup_chromo, num_of_bits);
            female_mask = ones(num_coup_chromo, num_of_bits);            
            
            %%%%%%%%%%%% BUG HERE - NEED TO FIX %%%%%%%%%%%%%%%
            male_mask(:, 1:crossover_point) = 0;
            female_mask(:,crossover_point + 1: num_of_bits) = 0;
            tmp_male = male .* male_mask;
            tmp_female = female .* female_mask;
            new_male = tmp_male .+ tmp_female;
            
            % Compute for new female
            male_mask = ones(num_coup_chromo, num_of_bits);
            female_mask = ones(num_coup_chromo, num_of_bits);            
            male_mask(:, crossover_point + 1:num_of_bits) = 0;
            female_mask(:, 1:crossover_point) = 0;
            tmp_male = male .* male_mask;
            tmp_female = female .* female_mask;
            new_female = tmp_male .+ tmp_female;
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
            [obj.init_pop_left_dec_lst, obj.init_pop_right_dec_lst] = convert_bin_2_dec(obj, obj.init_pop_left_bin_lst, obj.init_pop_right_bin_lst);
                            
            % TO-DO: Compute the final F(x, y) values
            % Insert F function here.
            % f_value_init_lst will be update after each iteration
            obj.f_value_init_lst = compute_func(obj, obj.init_pop_left_dec_lst, obj.init_pop_right_dec_lst, obj.resolution); 
            obj = crossover_good_couples(obj);
            disp('Final f_value: ');
            disp(obj.f_value_init_lst);
        end
        
        
        %%%% RANDOM SELECT COUPLES %%%%%%%%%%%%%%
        function obj = crossover_good_couples(obj)
            %%% TO-DO: Need to calculate based on cross probability
            %obj_use_to_cross always need to even.
            num_of_chromosome_4_crossover = obj.cross_prob * obj.pop_size / 100
            
            for iter = 1:1:obj.max_iter
                printf('Iteration #%d\n', iter);
                %%% Step 1: Initialize and select random couples base on the crossover probability.
                rand_idx = randperm(obj.pop_size, num_of_chromosome_4_crossover);
                crossover_points = randperm(obj.left_digits_num + obj.right_digits_num - 1, num_of_chromosome_4_crossover * 0.5);
                disp(crossover_points);
                couples_selected_lst = zeros(1, num_of_chromosome_4_crossover);
                chromosome_bin_selected = zeros(num_of_chromosome_4_crossover, obj.left_digits_num + obj.right_digits_num);
                
                couples_selected_lst = obj.f_value_init_lst(rand_idx, 1);
                
                % initialize chromosome
                init_bin_chromosome = [obj.init_pop_left_bin_lst, obj.init_pop_right_bin_lst];
                chromosome_bin_selected = init_bin_chromosome(rand_idx, :);
                
                %%% Step 2: Crossover couples
                male = chromosome_bin_selected(1:num_of_chromosome_4_crossover*0.5, :);
                female = chromosome_bin_selected(1 + num_of_chromosome_4_crossover*0.5 : num_of_chromosome_4_crossover, :);        
                [new_male, new_female] = crossover_a_couple(obj, male, female, crossover_points);
             
                %%% Step 3: Divide the chromosome into 2 components to recompute the F value
                new_bin_chromosome = [new_male; new_female];
                new_bin_left = new_bin_chromosome(:, 1:obj.left_digits_num);
                new_bin_right = new_bin_chromosome(:, obj.left_digits_num + 1: obj.left_digits_num + obj.right_digits_num);
     
                % ---> Recompute F value
                [new_dec_left ,new_dec_right] = convert_bin_2_dec(obj, new_bin_left, new_bin_right);
                new_f_values = compute_func(obj, new_dec_left, new_dec_right, obj.resolution);
                obj.new_f_value_lst = obj.f_value_init_lst;
                
                %%% Step 4: Compare the new F value with old F value before crossover, if bigger keep new and vice versa. 
                for i = 1:1:num_of_chromosome_4_crossover
                    if new_f_values(i, 1) > obj.new_f_value_lst(rand_idx(i), 1)
                        obj.new_f_value_lst(rand_idx(i),1) = new_f_values(i, 1);
                        
                        %Push back data to init data to store newest valid pop
                        obj.init_pop_left_bin_lst(rand_idx(i),:) = new_bin_left(i, :);
                        obj.init_pop_right_bin_lst(rand_idx(i), :) = new_bin_right(i, :);
                        obj.init_pop_left_dec_lst(rand_idx(i), 1) = new_dec_left(i, 1);
                        obj.init_pop_right_dec_lst(rand_idx(i), 1) = new_dec_right(i, 1);
                    end
                end
                obj.f_value_init_lst = obj.new_f_value_lst;
            end
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