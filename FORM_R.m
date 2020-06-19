function [u_standard, beta_final, p_f] = FORM_R(g, variables_list, dist, init_search_point)
%% Failure probability and the realibiliaty index using the HLRF for implicit limit state functions
%{
MIT License

Copyright (c) 2020 Iago Pereira Lemos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
---------------------------------------------------------------------------
Created by: Iago Pereira Lemos (lemosiago123@gmail.com)
Federal University of Uberlândia, School of mechanical engineering
---------------------------------------------------------------------------
Inputs
g: limit state function
variables_list: symbolic array with the variables in the limit state
function
dist: array with the distribution objects of which variable in
variables_list, following its order
init_search_point: array with the first values for each variable to start
searching for the solution

Outputs
u_standard: Coordinates of the design point in the standard space
beta_final: Reliability index value
p_f: Probability of failure
%}
%% Defining the gradient of the limit state equation
dg = gradient(g, variables_list);

%% Initializations 
maxit = 50;                               %maximum number of iterations    
beta  = ones(1, maxit);                   %initializing beta array for comparison  
beta(1) = 3;                              %first value of Beta
search_point = init_search_point;         %taking the initial search point into the current search point
cosines = ones(1, length(variables_list));%first values for the directions cosines
%% Solver
k = 1;
it = 0;
while true
    it = it + 1;
    % 1 step: Verify which variable is not normal and transform it to
    %considering the search point using the Rosenblatt transformation
    j = 0;
    for i=1:1:length(dist)
        if strcmp(class(dist(i)), class(makedist('Normal'))) == 0
            j = j + 1;
            non_normal2normal(j) = Rosenblatt_transform(dist(i), search_point(i));
        end
    end
    
    %initializing the normal means and stds vectors
    N_means = ones(1, length(dist));  
    N_std   = ones(1, length(dist));
    kk = 0;
    
    for i=1:1:length(dist)
        if strcmp(class(dist(i)), class(makedist('Normal'))) == 1
            N_means(i) = mean(dist(i));
            N_std(i)   = std(dist(i));
        elseif strcmp(class(dist(i)), class(makedist('Normal'))) == 0
            kk = kk + 1;
            N_means(i) = mean(non_normal2normal(kk));
            N_std(i)   = std(non_normal2normal(kk));
        end
    end
   
    %2 step: Computing the derivate in the current search point 
    derivate = double(subs(dg, variables_list, search_point))';
    
    %3 step: Computing direction cosines at the design point
    new_cosines  =(derivate.*N_std)/(sqrt(sum((derivate.*N_std).^2)));
    
    %4 step: Computing new search point
    new_search_point = N_means - (new_cosines.*N_std)*beta(k);
    search_point = new_search_point;
    
    %5 step: Verifying if the cosines converges
    if max(abs(new_cosines - cosines)) <= 0.005
        k = k+1;
        syms b %b = beta
        %6 step: Computing new design point with unknown beta
        %and solving g for this new design point
        design_point = N_means - b*(new_cosines.*N_std);
        g_value = subs(g, variables_list, design_point);
        beta_val = double(solve(g_value == 0, b));
        %7 step: Updating Beta
        beta(k) = double(min(beta_val));
        %Verifying if the new beta converges
        if abs(beta(k) - beta(k-1)) <= 0.00001
            break;
        end
    end
    cosines = new_cosines;
end

%% Results
beta_final = beta(k);
design_point_final = N_means - beta_final*(new_cosines.*N_std);
u_standard = N_means + design_point_final.*N_std;   %computing the final design point in the standard space
p_f = normcdf(-beta_final);                         %computing the failure probability
n_iterations = it;
%% Printing results
fprintf('Using the Rackwitz algorithm for implicit limit state functions:\n');
fprintf('Iterations: %g\nReliability index = %g\nFailure probability = %g\n\n', n_iterations, beta_final, p_f);

return
