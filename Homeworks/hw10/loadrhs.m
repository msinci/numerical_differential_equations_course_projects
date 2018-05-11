function b = loadrhs(fun,m,Lx,n,Ly)
%
% Creates the right-hand side vector for a linear PDE problem 
%
%               (PDE operator)u = f
%
% with Dirichlet boundary conditions on an nxm grid on a 
% rectangular domain [0,Lx]x[0,Ly]. 
%
% The inputs are 
% Required:
%   fun     = user-created function handle (@functionname) for 
%             the right-hand side function f; this should return 
%             f(x,y) for scalar inputs x and y. 
%   m       = number of meshpoints in the x-direction. 
%
% Optional: 
%   Lx      = length of rectangle in the x-direction (default 1). 
%   n       = number of meshpoints in the y-direction (default m).
%   Ly      = length of rectangle in the y-direction (default 1 or Lx). 
%
% The output is 
%   b       = right-hand side vector 
% 

if nargin < 2, 
  error('loadrhs: The number of arguments must be => 2; see comments.'); 
elseif nargin > 5, 
  error('loadrhs: The number of arguments must be <= 5; see comments.'); 
end 
if nargin == 2
  Lx = 1; n = m; Ly = 1;
elseif nargin == 3 
  n = m; Ly = Lx;
elseif nargin == 4
  Ly = Lx;
end 
if isempty(Lx), Lx = 1; end
if isempty(n), n = m; end
if isempty(Ly), Ly = Lx; end

hx = (m+1)\Lx;
hy = (n+1)\Ly; 

x = hx:hx:1-hx;
y = hy:hy:1-hy;

b = zeros(m,n);
for j = 1:n
    for i = 1:m
      b(i,j) = feval(fun,x(i),y(j));
    end
end
b = reshape(b,m*n,1);
