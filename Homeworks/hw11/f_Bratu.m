function F_val = f_Bratu(u)
%
% Evaluates F(u) for the Bratu problem 
%
%   (Laplacian)u + lambda*exp*u = 0 in D = [0,1]x[0,1] 
%                             u = 0 on the boundary of D
% 
% The discretization is by finite differences on an mxm regular 
% grid of interior points in [0,1]x[0,1]. 
%

global L_GL; 
global lambda_GL;

% lambda_GL is passed in from newton_bratu_demo.

n = size(u,1);
m = sqrt(n); 
if n ~= fix(m)^2, 
  error('n must be m^2, where m is the number of gridpoints per side.'); 
end

if size(L_GL,1) ~=n, 
  L_GL = laplacian(m);
end

F_val = L_GL*u + lambda_GL*exp(u);

