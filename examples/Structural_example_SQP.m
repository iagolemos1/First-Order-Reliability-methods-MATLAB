clear, close all, clc
%% Structural example: The moment capacity of a singly reinforced rectangular prism concrete beam - Normal variables
%{ 
Original equation:
As*fy*d*(1 - (neta*As*fy/(b*d*fc)) = M
g() = As*fy*d*(1 - (neta*As*fy/(b*d*fc)) - M

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

%% Writing g() in the explicit form (standard space)
g = @(x) (x(1,:)*0.05616 + 1.56)*(x(2,:)*7.155 + 47.7)*(x(3,:)*1.1352 + 13.2)*(1 - ((x(4,:)*0.0295 + 0.59)*(x(1,:)*0.05616 + 1.56)*(x(2,:)*7.155 + 47.7)/((x(5,:)*0.36 + 8)*(x(3,:)*1.1352 + 13.2)*(x(6,:)*0.735 + 3.5)))) - (x(7,:)*55.4625 + 326.25);

%% Solving
FORM_SQP(g,7);