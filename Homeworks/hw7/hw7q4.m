% Mehmet Sinan INCI 
% MA 512 - HW 7 Q4
clc;
clear all;  % clear variables
close all;  % close old figures
format long;

%% Now illustrate setting tolerances using ode15s options
% initial values
y_10 = 1;
y_20 = 1.8;
y_30 = 1.8;
tinit = 0;
tfinal = 290;

% Set the options for ode solvers
opt1 = odeset('Stats','on'); % show # of steps etc
% ode15
[tvals, yvals] = ode15s('q4f',[tinit tfinal],[y_10 y_20 y_30],opt1);
nsteps = size(tvals,1); % calculate the n
fprintf('nn nsteps = %dnn', nsteps);
figure(1);
plot3(log(yvals(:,1)),log(yvals(:,2)),log(yvals(:,3)));
grid on
