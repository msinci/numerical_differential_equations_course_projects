function [A,H,V] = condif(c,d,e,m,Lx,n,Ly)
%
% Returns the sparse matrix for the discrete convection-diffusion 
% operator 
%
%             -(Laplacian)u + cu + du_x + eu_y 
%
% with Dirichlet boundary conditions on an nxm grid on a 
% rectangular domain [0,Lx]x[0,Ly]. 
%
% Optionally returns the "horizontal" and "vertical" operators 
% for use with ADI. 
%
% The inputs are as follows: 
%
% Required:
%   c         = coefficient of zero-order term. 
%   d         = coefficient of x-derivative term. 
%   e         = coefficient of y-derivative term. 
%   m         = number of meshpoints in the x-direction. 
%
% Optional: 
%   Lx        = length of rectangle in the x-direction (default 1). 
%   n         = number of meshpoints in the y-direction (default m).
%   Ly        = length of rectangle in the y-direction (default 1 or Lx). 
%
% The outputs are  
%   A         = sparse discrete convection-diffusion operator. 
%   H	      = "horizontal" operator.
%   V	      = "vertical" operator.
% 
if nargin < 4, error('Number of arguments must be => 4.'); end 
if nargin == 4,
  n = m; Lx = 1; Ly = 1;
end
if nargin == 5, 
  n = m; Ly = Lx;
end
if nargin == 6,
  Ly = Lx;
end
if nargin > 7, error('Number of arguments must be <= 4.'); end 

hx = (m+1)\Lx;
hy = (n+1)\Ly; 

Im = speye(m,m);
Em = sparse(2:m,1:m-1,1,m,m);
DDm = -2*Im + Em + Em';

In = speye(n,n);
En = sparse(2:n,1:n-1,1,n,n);
DDn = -2*In + En + En';

Sm = ((hx^2)\1)*DDm + ((2*hx)\d)*(Em' - Em) + (2\c)*Im;
Sn = ((hy^2)\1)*DDn + ((2*hy)\e)*(En' - En) + (2\c)*In;

H = kron(In,Sm);
V = kron(Sn,Im); 

A = V + H; 
