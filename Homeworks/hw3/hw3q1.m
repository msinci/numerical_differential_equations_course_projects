% Mehmet Sinan INCI 
% MA 512 - HW 3 Q1

%% part 1
clc;
clear all;  % clear variables
close all;  % close old figures

tspan = [0 20]; % our given range
% y = ones(20,3);
y_0 = [1;1];    % initial values
options = odeset('RelTol',10^-8, 'AbsTol', 10^-8); % setting paremeters for ode solver
[t, y] = ode45(@q1f, tspan,y_0,options); % using the ode45 solver

step_nr = size(t,1); % number of steps, we will print this value
fprintf('\n step_nr = %d, no. f-evals = %d',step_nr, 6*step_nr);

figure; plot(t,y);
figure; plot(y(:,1),y(:,2));
axis([-1.5 1.5 -1.5 1.5]); axis square;
