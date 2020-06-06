clear, close all, clc
%% Structural example: The moment capacity of a singly reinforced rectangular prism concrete beam - Normal variables
%{ 
Original equation:
As*fy*d*(1 - (neta*As*fy/(b*d*fc)) = M
g() = As*fy*d*(1 - (neta*As*fy/(b*d*fc))) - M

As: area of the tension reinforcing bars
fy: yield stress of the reinforcing bars
d : distance from the extreme compression fiber to the centroid of the
tension reinforcing bars
neta : the concrete stress block parameter
fc   : the compressive strength of concrete
b    : of the compression face of the member.
M    : the moment capacity or resistance

Variables values table

Rand. Variables    Mean     Standard deviation
   As (in^2)       1.56           0.05616
   fy (ksi)        47.7           7.155
   d  (in)         13.2           1.1352
   neta            0.59           0.0295
   fc (ksi)        3.5            0.735
   b  (in)         8.0            0.36
   M  (kip-in)     326.25         55.4625
%}  

% Writing g function
syms As fy d neta fc b M
g = As*fy*d*(1 - (neta*As*fy/(b*d*fc))) - M;

%% 1 - Considering all normal variables

variables_list1 = [As, fy, d, neta, fc, b M];  %list of varibles in the right order

% Creating dists objects
As_dist1   = makedist('Normal', 'mu', 1.56,'sigma', 0.05616);
fy_dist1   = makedist('Normal', 'mu', 47.7,'sigma', 7.155);
d_dist1    = makedist('Normal', 'mu', 13.2,'sigma', 1.1352);
neta_dist1 = makedist('Normal', 'mu', 0.59,'sigma', 0.0295);
fc_dist1   = makedist('Normal', 'mu', 3.5,'sigma', 0.735);
b_dist1    = makedist('Normal', 'mu', 8.0,'sigma', 0.36);
M_dist1    = makedist('Normal', 'mu', 326.25,'sigma', 55.4625);

% Creating the vector of distribution objects
dist1 = [As_dist1, fy_dist1, d_dist1, neta_dist1, fc_dist1, b_dist1, M_dist1];

% Initial search points
init_search_point1 = [1,1,1,1,1,1,1];


%Printing results
fprintf('------Considering all normal variables:------\n')
FORM_R(g, variables_list1, dist1, init_search_point1);
FORM_RF(g, variables_list1, dist1, init_search_point1);
fprintf('\n')

%% 2 - Considering M as lognormal variable
%Obs: It's easier considering M as lognormal because it's possible to
%get the lognormal distribution by knowing the mean and the std of the
%variable.

M_mean = 326.25;       %M variable mean
M_std  = 55.4625;      %M variable std
cov_M  = M_std/M_mean; %Covariation of the variable M

M_log_std = sqrt(log(1 + cov_M^2));           %Computing the std of the lognormal dist
M_log_mean = log(M_mean) - 0.5*(M_log_std^2); %Computing the mean of the lognormal dist

%Creating lognormal dist
M_dist2    = makedist('Lognormal', 'mu', M_log_mean,'sigma', M_log_std);

variables_list2 = [As, fy, d, neta, fc, b M];

dist2 = [As_dist1, fy_dist1, d_dist1, neta_dist1, fc_dist1, b_dist1, M_dist2];

init_search_point2 = [1,1,1,1,1,1,1];

%Printing results
fprintf('------Considering M as a lognormal variable:------\n')
FORM_R(g, variables_list2, dist2, init_search_point2);
FORM_RF(g, variables_list2, dist2, init_search_point2);
fprintf('\n')