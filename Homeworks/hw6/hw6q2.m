% Mehmet Sinan INCI 
% MA 512 - HW 6 Q2

%% 
clc;
clear all;  % clear variables
close all;  % close old figures

tspan = [0 20]; % our given range
% y = ones(20,3);
y_0 = [1;1];    % initial values

% options = odeset('Jacobian',@J,'RelTol',10^-8, 'AbsTol', 10^-8, 'Stats', 'on'); % setting paremeters for ode solver

options = odeset('RelTol',10^-8, 'AbsTol', 10^-8, 'Stats', 'on'); % setting paremeters for ode solver
[t, y] = ode113(@q2f, tspan,y_0,options); % using the ode45 solver

% figure; plot(t,y);
% figure; plot(y(:,1),y(:,2));
% axis([-1.5 1.5 -1.5 1.5]); axis square;
