function fval = rhsfun(x,y)
%
% Evaluates a function f(x,y) for specified x and y; called by loadrhs.m. 
%
% The inputs are 
%   x       = scalar x argument. 
%   y       = scalar y argument. 
%
fval = 2*(x*(1-x) + y*(1-y)); 
