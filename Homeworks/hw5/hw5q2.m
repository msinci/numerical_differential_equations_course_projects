% Mehmet Sinan INCI 
% MA 512 - HW 5 Q2

clc;
clear all;  % clear variables
close all;  % close old figures
format long; 
Low_bound = 0;  % lower bound for function
Up_bound = pi;  % upper bound for function
index=1;    % start index for the vector

for n = [5, 50, 500]
    y(1) = 0; % initial value
    h = pi/n; % step size = upper bound-lower bound / n
    tra = [0:h:pi]; % area of operation
    exact = sin(tra);
    for i=1:n % try all 3 number of steps
        t = tra(i); % holding the t value
        tp = tra(i+1); % holding the next t value
        part1 = (1-500*h) * y(i)/(1+500*h);
        part2 = h/(2+1000*h)*(1000*sin(tra(i))+cos(tra(i))+1000*sin(tra(i+1))+cos(tra(i+1)));
        y(i+1) = part1 + part2;
    end
    err(index) = max(abs(y-exact)); % calculate the error
    index = index+1;    % increment the index in each iteration
end

err
