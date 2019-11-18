function output = genetic_optimization(low1, high1, low2, high2, pop_size, resolution, max_iter)
    if nargin == 4 
        warning('[WARNING] Not determine the pop size, resolution, max_iter');
        fprintf('The default pop size will be: 100\n');
        fprintf('The default resolution will be: 3\n');     
        fprintf('The default max iteration will be: 100\n');   
        pop_size = 100;
        resolution = 3;
        max_iter = 100;
    elseif nargin == 5
        warning('[WARNING] Not determine the resolution');
        fprintf('The default resolution will be: 3\n');
        fprintf('The default max iteration will be: 100\n');
        resolution = 3;
        max_iter = 100;
    elseif nargin == 6
        warning('[WARNING] Not determine the max iteration');
        fprintf('The default max iteration will be: 100\n');
        max_iter = 100;
    elseif nargin < 4
        fprintf('[ERROR]Not enough input arguments\n');
        return;
    elseif nargin > 7
        fprintf('Redundance input arguments\n');
        return;
    end
    
    [parent_pop, left_elems, right_elems] = initialize_pop(low1, high1, low2, high2, pop_size, resolution, max_iter);
    natural_selection(parent_pop, left_elems, right_elems, pop_size);
end
