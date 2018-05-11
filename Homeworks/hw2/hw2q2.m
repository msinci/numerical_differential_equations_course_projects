% Mehmet Sinan INCI 
% MA 512 - HW 2 Q2

%%
clc;
clear all;
close all;
N = 100;  % we will try different step sizes to check error
h=pi/N;            % step size
t = 0:h:pi;
y = zeros(1,length(t)); 
y(1) = 0;          % initial condition
exact(1) = 0;
error(1) = 0;

F_xy = @(t,y) exp(t).*cos(t)+y; % approx. function

for i=1:(length(t)-1)                      % calculation loop
    V_1 = F_xy(t(i),y(i));                      % calculate V1
    V_2 = F_xy(t(i)+0.5*h,y(i)+0.5*h*V_1);      % calculate V2
    V_3 = F_xy((t(i)+0.5*h),(y(i)+0.5*h*V_2));  % calculate V3
    V_4 = F_xy((t(i)+h),(y(i)+V_3*h));          % calculate V4
    y(i+1) = y(i) + (1/6)*(V_1+2*V_2+2*V_3+V_4)*h;  % calculate the final result
end

exact = exp(t).*sin(t); % exact solution
error = abs(exact - y); % error

N      % printing out the N value
h  % printing out the step size
max_error = max(error) % printing out the max error
% 
% plot(y);title('approximation');
% hold on;
% plot(exact); hold off;
% figure; plot(error);title('error');
