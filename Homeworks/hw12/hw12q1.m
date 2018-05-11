% Mehmet Sinan INCI 
% MA 512 - HW 12 Q1
clc;
clear all;  % clear variables
close all;  % close old figures
format long; 

tint = 0:10;
lamda = 3.51;   %% we try lamda=3.51 and 3.52

y_10 = 1;
y_20 = 1.8;
y_30 = 1.8;
tinit = 0;
tfinal = 10;
option = odeset('Stats','on'); % show # of steps etc
[tvals, yvals] = ode15s('mol_demo_1D',[tinit tfinal],[y_10 y_20 y_30],option);
nsteps = size(tvals,1); % calculate the n
fprintf('nn nsteps = %dnn', nsteps);
figure(1);
plot3(log(yvals(:,1)),log(yvals(:,2)),log(yvals(:,3)));
grid on


%%
% plot(T,max(abs(U'))');