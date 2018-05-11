function [x,itno,resvec] = gaussseidel(x, A, b, tol, itmax)
%
% Implements the Gauss-Seidel method to produce x such that ||b-Ax|| <= tol. 
%
% The inputs are 
%   x       = initial iterate 
%   A       = coefficient matrix
%   b       = right-hand side
%   tol     = stopping tolerance
%   itmax   = maximum allowable number of iterations 
%
% The returns are 
%   x       = final iterate
%   itno    = vector of iteration numbers
%   resvec  = vector of residual norms
%
L = tril(A);
U = triu(A,1);
itno = [0];
rsnrm = norm(b - A*x);
resvec = [rsnrm];
j = 0;
while rsnrm > tol
   j = j + 1;
   if j > itmax, break, end 
   x = L\(b - U*x);
   itno = [itno,j];
   rsnrm = norm(b - A*x); 
   resvec = [resvec,rsnrm];
end
