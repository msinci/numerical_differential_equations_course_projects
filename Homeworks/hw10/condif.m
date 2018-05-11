function [A,H,V] = condif(m,coeffs,Lx,n,Ly)
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
% H and V for use with ADI. 
%
% The inputs are as follows: 
%
% Required:
%   m         = number of meshpoints in the x-direction. 
%   coeffs    = row vector of coefficients [c d e].
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

if nargin < 2, 
  error('condif: Number of arguments must be => 2.'); 
elseif nargin > 5, 
  error('condif: Number of arguments must be <= 5.'); 
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

c = coeffs(1,1); d = coeffs(1,2); e = coeffs(1,3); 

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
