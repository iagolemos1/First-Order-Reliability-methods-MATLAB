function [u_standard, beta_final, p_f] = FORM_RF(g, variables_list, dist, init_search_point)
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
search_point = init_search_point; %taking the first point of searching as the user's initial search point
maxit = 50;                       %maximum number of iterations
j = 1;                            %counter for iterations
beta = ones(1,maxit);             %initializing beta array for comparison

%% Solver

while true
    % 1 step: Verify which variable is not normal and transform it to
    %considering the search point using the Rosenblatt transformation
    k = 0;
    for i=1:1:length(dist)
        if strcmp(class(dist(i)), class(makedist('Normal'))) == 0
            k = k + 1;
            non_normal2normal(k) = Rosenblatt_transform(dist(i), search_point(i));
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
    
    %2 step: Computing the new design point in the standard space
    new_design_point_standard = (search_point - N_means)./N_std;
    design_point_standard = new_design_point_standard;
    
    %3 step: Computing the derivate in the current search point 
    derivate = double(subs(dg, variables_list, search_point));
    partial_derivate = derivate'.*N_std;
    
    %4 step´: Computing the new design point in the standard space
    esc1 = 1/sum(partial_derivate.*partial_derivate);
    esc2 = sum(partial_derivate.*design_point_standard) - double(subs(g, variables_list, search_point));
    esc3 = esc1 * esc2;
    new_design_point_standard = esc3*partial_derivate;
    
    %5 step: Computing new beta
    beta(j+1) = norm(new_design_point_standard);
    
    %6 step: Computing the new search point (in the original space)
    new_search_point = new_design_point_standard.*N_std + N_means;
    
    %7 step: Computing the g() in the new search point
    g1 = double(subs(g, variables_list, new_search_point));
    
    search_point = double(new_search_point);
    
    %8 step: Verifying convergence
    if (abs(beta(j+1) - beta(j))<=0.001) && (abs(g1)<=0.001)
        break;
        
    else
        j=j+1;
    end
end

%% Results
beta_final   = beta(j+1);                 %taking the final beta
u_standard   = new_design_point_standard; %taking the new point in the standard space
n_iterations = j;                         %taking the number of iterations    
p_f = normcdf(-beta_final);               %taking the failure probability

%% Printing results
fprintf('Using the Rackwitz-Fiessler algorithm for implicit limit state functions:\n');
fprintf('Iterations: %g\nReliability index = %g\nFailure probability = %g\n\n', n_iterations, beta_final, p_f);

return
