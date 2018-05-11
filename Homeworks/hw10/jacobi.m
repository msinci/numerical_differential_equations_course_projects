function [x,itno,resvec] = jacobi(x, A, b, tol, itmax)
%
% Implements Jacobi iteration to produce x such that ||b-Ax|| <= tol. 
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
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
C = L+U; 
itno = [0];
rsnrm = norm(b - A*x);
resvec = rsnrm;
j = 0;
while rsnrm > tol
   j = j + 1;
   if j > itmax, break, end 
   x = D\(b - C*x);
   itno = [itno,j];
   rsnrm = norm(b - A*x); 
   resvec = [resvec,rsnrm];
end
