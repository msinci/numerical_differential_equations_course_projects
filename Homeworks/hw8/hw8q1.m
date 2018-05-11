% Mehmet Sinan INCI 
% MA 512 - HW 7 Q1
clc;
clear all;  % clear variables
close all;  % close old figures
format long; 

% tspan = [0 100]; % our given range
h = [0.1, 0.01, 0.001];

Low_bound = 0;  % lower bound for function
Up_bound = 100;  % upper bound for function
A = [0,1;-1,0];
y(:,:,1) = [1,0]';    % initial values

% The midpoint method

for i = 1:length(h) % loop runs 3 times for each step size
    t = (Low_bound:h(i):Up_bound);
    % We first obtain the solution at x = h;
    y(:,:,2) = y(:,:,1) + h(i)*A*y(:,:,1);
    % Midpoint method
    for j = 2:length(t)-1
        y(:,:,j+1) = y(:,:,j-1) + 2*h(i)*A*y(:,:,j);
    end
    % Calculating the exact solution
    for k = 1: length(t)
        exact(:,:,k) = q1f( t(k) ); % use the exact solution function 
                                   % to calculate the exact value at points
    end
    err(i) = max(max(y - exact)); % We only record the maximum error 
                                 % and discard smaller errors
    fprintf('H is %d.nn, Max Error is %d\n',h(i), err(i));
end