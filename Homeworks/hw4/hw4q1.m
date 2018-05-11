% Mehmet Sinan INCI 
% MA 512 - HW 4 Q1

clc;
clear all;  % clear variables
close all;  % close old figures
tspan = [0 20];
yvals = zeros(20,2);

yvals(1,1) = 0; % Q(t)
yvals(1,2) = 2; % Q'(t)
% [tvals, yvals] = ode45(@q1f, tspan,[0;2]); % using the ode45 solver
options = odeset('RelTol',10^-12, 'AbsTol', 10^-12); % setting paremeters for ode solver
[tvals, yvals] = ode45(@q1f, tspan,[0;1.999],options); % with options
plot(tvals, yvals(:,1));
hold on;
yvals = zeros(20,2);
[tvals, yvals] = ode45(@q1f, tspan,[0;2],options); % with options
plot(tvals, yvals(:,1));
hold on;
yvals = zeros(20,2);
[tvals, yvals] = ode45(@q1f, tspan,[0;2.001],options); % with options
plot(tvals, yvals(:,1));

