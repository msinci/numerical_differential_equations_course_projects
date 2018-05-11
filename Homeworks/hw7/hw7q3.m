% Mehmet Sinan INCI 
% MA 512 - HW 7 Q3
clc;
clear all;  % clear variables
close all;  % close old figures
format long;

y_10 = 1; % x initial value
y_20 = 1; % x'Initial value
t_0 = 0;
t_f = 20;

% Since we are running multiple ODE solvers with same options, 
% it makes sense to save options and feed it into these solvers
opt1 = odeset('RelTol',1e-8,'AbsTol',1e-8,'Jacobian',[0 1;-1 -100],'Stats','on');
opt2 = odeset('RelTol',1e-8,'AbsTol',1e-8,'Jacobian',[0 1;-1 -100],'BDF','on','Stats','on');
opt3 = odeset('RelTol',1e-8,'AbsTol',1e-8,'Stats','on'); % ode45

fprintf('\node45 solver results\n');
[tvals1,yvals1] = ode45(@q3f,[t_0 t_f],[y_10 y_20],opt3);

fprintf('\node113 solver results\n');
[tvals2,yvals2] = ode113(@q3f,[t_0 t_f],[y_10 y_20],opt1);

% Using BDF
fprintf('\node15 solver results using BDF\n');
[tvals3,yvals3] = ode15s(@q3f,[t_0 t_f],[y_10 y_20],opt2);
% Using NDF
fprintf('\node15 solver results using NDF\n');
[tvals4,yvals4] = ode15s(@q3f,[t_0 t_f],[y_10 y_20],opt1);
