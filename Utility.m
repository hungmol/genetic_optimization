classdef Utility
    properties (SetAccess = public)
        test = 10;
    end
    methods(Static)
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
        
    end
  
endclassdef
