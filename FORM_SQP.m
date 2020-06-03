function [u_standard, beta, p_f] = FORM_SQP(g, d)
%% Failure probability and the realibiliaty index using the Sequential Quadratic Programming solver
%{
---------------------------------------------------------------------------
Created by: Iago Pereira Lemos (lemosiago123@gmail.com)
Federal University of Uberlândia, School of mechanical engineering
---------------------------------------------------------------------------
Inputs
g: limit state function
d: dimensions of the limit state function

Outputs
u_standard: Coordinates of the design point in the standard space
beta: Reliability index value
p_f: Probability of failure
%}

%% Defining the objective function
ob_func = @(u) norm(u); 

%% Defining the parameters of the fmincon function
u_initial  = repmat(0.1, d, 1); %returns a matrix of 0.1 in the for d x 1,
%defining the initial search point

% The constraints are an empty array because we use a non linear function
% in order to define the constraints.
A          = [];                % linear equality constraints
b          = [];                % linear equality constraints
Aeq        = [];                % linear inequality constraints
beq        = [];                % linear inequality constraints
lb         = [];                % lower bound constraints
ub         = [];                % upper bound constraints

% Defining the constraints: cons_func(u) = 0
cons_func = @(u) deal([], g(u)); 

%% Defining solver and geting the solution

% Defining solver options, note: using sequential quadratic programming
%'sqp'. 
solver_options = optimoptions('fmincon','Display','off','Algorithm','sqp');

% Using fmincom with the whole defined parameters
[u_standard,beta,~,output] = fmincon(ob_func, u_initial, A, b, Aeq, beq, lb, ub, cons_func, solver_options);
n_iterations = output.iterations; %geting the number of interations
algorithm  = output.algorithm; %geting the algorithm 

% computing the failure probability
p_f = normcdf(-beta);

% printing results
fprintf('Using fmincon function with the %s solver:\n',algorithm);
fprintf('Iterations: %g\nReliability index = %g\nFailure probability = %g\n\n', n_iterations, beta, p_f);

return;

