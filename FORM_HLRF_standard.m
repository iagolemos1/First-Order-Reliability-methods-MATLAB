function [u_standard, beta, p_f] = FORM_HLRF_standard(g,variables_list)
%% Failure probability and the realibiliaty index using the HLRF method for standard limit state functions 
%{
---------------------------------------------------------------------------
Created by: Iago Pereira Lemos (lemosiago123@gmail.com)
Federal University of Uberlândia, School of mechanical engineering
---------------------------------------------------------------------------
Inputs
g: limit state function
variables_list: vector of symbolic variables in the limit state function

Outputs
u_standard: Coordinates of the design point in the standard space
beta: Reliability index value
p_f: Probability of failure
%}

%% initialization
tol     = 0.001;                %error
maxiter = 50;                   %maximum number of iterations
d       = length(variables_list)%number of variables in the l.s.f.
x       = ones(d,maxiter);      %initialization of the vector of design points
%% HLRF method
dg = gradient(g, variables_list);
k = 1;  %initializing the counter
while true
 
   dg_x_k      = subs(dg, variables_list, x(:,k));     %comuting the gradient of the l.s.f. in the k design point
   g_x_k       = subs(g, variables_list, x(:,k));      %computign the l.s.f value in the k design point
   norm_dg_x_k = norm(dg_x_k);                         %computing the norm of the grandient of l.s.f in the k design point
   
   
   % doing some maths
   esc_1      = 1/norm_dg_x_k^2; 
   esc_2      = (dg_x_k'*x(:,k) - g_x_k);
   esc_3      = esc_1 * esc_2;
   
   x(:,k+1)   = esc_3 * dg_x_k;
   
   % next iteration
   if (norm(x(:,k+1)-x(:,k)) <= tol)  || (norm_dg_x_k <= tol)
      break;
   else
      k = k+1;
   end
end


beta = norm(x(:,k));

n_iterations = k;
p_f = normcdf(-beta);
u_standard = x(:,k);

fprintf('Using the Hasofer-Lind-Rackwitz-Fiessler algorithm for standard limit state functions:\n');
fprintf('Iterations: %g\nReliability index = %g\nFailure probability = %g\n\n', n_iterations, beta, p_f);

return;