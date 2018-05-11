function [yvals,tvals] = forward_euler(Fname,y_0,t_0,t_f,nsteps)
%
% This function implements the forward Euler method to solve 
% the initial value problem 
%
%       y' = f(t,y),   y(t_0) = y_0
%
% over the interval [t_0,t_f] using a constant step-size. This 
% can be either a scalar problem or a vector system. 
%
% The inputs are 
%  Fname   = function handle (@fname) specifying an M-file for 
%	     evaluating f(t,y). The M-file should begin with 
%	     "function fval = fname(t,y)," where fval is the 
%	     value of f(t,y) computed by the file. 
%  y_0     = initial value of y(t).
%  t_0     = initial t. 
%  t_f     = final t. 
%  nsteps  = number of steps across [t_0,t_f].
%
% The output is 
%  yvals   = array of dimension (length of y)x(nstep+1) whose 
%	     columns are the computed approximations of y at 
%	     the time values. 
% tvals    = vector of time values. 
%
% Note: y_0, y, and f-values returned by @fname are expected to 
% be scalars or column vectors of the same length. 
%

h = nsteps\(t_f - t_0); 
tvals = [t_0:h:t_f]; 
yvals = zeros(size(y_0,1),nsteps + 1); 
yvals(:,1) = y_0; 
for i = 1:nsteps 
  yvals(:,i+1) = yvals(:,i) + h*feval(Fname,tvals(1,i),yvals(:,i));
end
