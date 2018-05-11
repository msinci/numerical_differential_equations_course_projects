% Mehmet Sinan INCI 
% MA 512 - HW 5 Q2

clc;
clear all;  % clear variables
close all;  % close old figures
tspan = [0 20];

options = odeset('RelTol',10^-8, 'AbsTol', 10^-8); % setting paremeters for ode solver
[tvals, yvals] = ode45(@q2f, tspan, [1.1;1;1],options); % with options
plot3(yvals(:,1), yvals(:,2), yvals(:,3));

