function [x,itno,resvec] = sor(x, A, b, w, tol, itmax)
%
% Implements the SOR method to produce x such that ||b-Ax|| <= tol.
%
% The inputs are 
%   x       = initial iterate 
%   A       = coefficient matrix, assumed to be sparse 
%   b       = right-hand side
%   w       = relaxation parameter
%   tol     = stopping tolerance
%   itmax   = maximum allowable number of iterations 
%
% The returns are 
%   x       = final iterate
%   itno    = vector of iteration numbers
%   resvec  = vector of residual norms
%
L = tril(A,-1);
U = triu(A,1);
D = diag(diag(A)); 
itno = [0];
rsnrm = norm(b - A*x);
resvec = rsnrm;
j = 0;
while rsnrm > tol
   j = j + 1;
   if j > itmax, break, end 
   x = (D + w*L)\(((1-w)*D - w*U)*x + w*b); 
   itno = [itno,j];
   rsnrm = norm(b - A*x); 
   resvec = [resvec,rsnrm];
end
