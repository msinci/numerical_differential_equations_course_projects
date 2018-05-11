% Mehmet Sinan INCI 
% MA 512 - HW 3 Q2

clc;
clear all;  % clear variables
close all;  % close old figures
tspan = [0 20];
options = odeset('RelTol',10^-6, 'AbsTol', 10^-12); % setting paremeters for ode solver
[t, y] = ode45(@hunt, tspan,[8;3],options); % using the ode45 solver

step = size(t,1) -1;
r = y(:,1);
f = y(:,2);

plot(r,f);


%% with event parameters

clc;
clear all;  % clear variables
close all;  % close old figures
tspan = [0 20];
options = odeset('Events',@hw3funevent);
[t,y,tEvent, yEvent] = ode45(@hunt, tspan,[8;3],options);
step = size(t,1)-1
r = y(:,1);
f = y(:,2);
plot(r,f);
