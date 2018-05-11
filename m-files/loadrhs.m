function b = loadrhs(fun,m,Lx,n,Ly) 
%
% Creates the right-hand side vector for a PDE problem Lu = f 
% on an mxn grid of interior points in a rectangular domain 
% [0,Lx]x[0,Ly]. 
%
% The inputs are as follows: 
%
% Required:
%   fun     = function handle (@functionname) for evaluating 
%	      the right-hand side function f. This should accept 
%	      scalar arguments x and y and return f(x,y). 
%   m       = number of meshpoints in the x-direction. 
%
% Optional: 
%   Lx      = length of rectangle in the x-direction (default 1). 
%   n       = number of meshpoints in the y-direction (default m).
%   Ly      = length of rectangle in the y-direction (default 1 or Lx). 
%
% The output is 
%   b       = right-hand side vector. 
% 
%
if nargin < 2, error('Number of arguments must be >= 2.'); end 
if nargin == 2,
  n = m; Lx = 1; Ly = 1;
end
if nargin == 3, 
  n = m; Ly = Lx;
end
if nargin == 4,
  Ly = Lx;
end
if nargin > 4, error('Number of arguments must be <= 4.'); end 

hx = (m+1)\Lx;
hy = (n+1)\Ly; 

x = hx:hx:1-hx;
y = hy:hy:1-hy;
for j = 1:n
    for i = 1:m
      b(i,j) = feval(fun,x(i),y(j));
    end
end
b = reshape(b,m*n,1);
