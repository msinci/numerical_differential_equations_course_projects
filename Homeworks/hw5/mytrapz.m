function [ out ] = mytrapz( func, a, b, n)
format long; 
% function takes the function f, lower bound a, upper bound b and n
    x = a; % initial x equal to lower bound, a
    h = pi/n;  % h value, the step size
    flow=func(a);   % calculate f at a, the lower bound
    for i = 1 : n-1
      x = x + h;    % increment by the step
      flow = flow + 2*func(x);  % calculate ned f(lower bound) value
    end
    flow = flow + func(b); % func(b) is fhigh
    out = (b - a) * flow/(2*n);     % return the integral value
end