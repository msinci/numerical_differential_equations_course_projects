function [L,H,V] = laplacian(m,Lx,n,Ly)
%
% Returns the sparse matrix for the discrete Laplacian operator 
% with Dirichlet boundary conditions on a regular nxm grid of
% interior points on a rectangular domain [0,Lx]x[0,Ly]. 
%
% Optionally returns the "horizontal" and "vertical" operators 
% H and V for use with ADI. 
%
% The inputs are as follows: 
%
% Required:
%   m            = number of meshpoints in the x-direction. 
%
% Optional: 
%   Lx          = length of rectangle in the x-direction (default 1). 
%   n           = number of meshpoints in the y-direction (default m).
%   Ly          = length of rectangle in the y-direction (default 1 or Lx). 
%
% The outputs are  
%   L           = sparse discrete Laplacian. 
%   H		= "horizontal" operator.
%   V		= "vertical" operator.
% 

if nargin > 4, error('laplacian: Number of arguments must be <= 4.'); end 

if nargin == 1,
  n = m; Lx = 1; Ly = 1;
elseif nargin == 2, 
  n = m; Ly = Lx;
elseif nargin == 3,
  Ly = Lx;
end
if isempty(Lx), Lx = 1; end
if isempty(n), n = m; end
if isempty(Ly), Ly = Lx; end

hx = (m+1)\Lx;
hy = (n+1)\Ly; 

Im = speye(m,m);
Em = sparse(2:m,1:m-1,1,m,m);
DDm = -2*Im + Em + Em';

In = speye(n,n);
En = sparse(2:n,1:n-1,1,n,n);
DDn = -2*In + En + En';

H = ((hx^2)\1)*kron(In,DDm);
V = ((hy^2)\1)*kron(DDn,Im);

L = H + V; 

