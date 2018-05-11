% Mehmet Sinan INCI 
% MA 512 - HW 2 Q1
%%
clc;
clear all;
close all;

% y' =y+(e^t *cos(t))
% initial condition
x(1) = 0;
y(1) = 0;
N = 100;
h = pi/N;% step size
t = 0:h:pi;

Fp = @(t,y) y+exp(t).*cos(t);

%Euler Forward
for i=1:(length(t)-1) 
    y(i+1) = y(i)+h*(y(i)+exp(t(i)).*cos(t(i)));
%     y(i+1) = y(i)+h*Fp(i,y(i));
    x(i+1) = x(i)+h;
end

exact = exp(t).*sin(t); % exact solution
error = abs(exact - y); % error
N      % printing out the N value
h      % printing out the step size
max_error = max(error) % printing out the max error

% plot
figure
set(gca,'Fontsize',15)
plot(x,y,'+-', 'Linewidth', 1)
xlabel('x')
ylabel('y')
