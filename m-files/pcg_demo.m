function pcg_demo(m,Pre_Name,Problem,fun,tol,itmax)
%
% Demonstrates the performance of the preconditioned conjugate 
% gradient method (PCG) on the problem 
%
%             -(Laplacian) u = f  or (Laplacian)^2 u = f  
%
% with Dirichlet boundary conditions on the unit square in R^2 
% discretized on an evenly spaced mxm mesh of interior points. 
%
% m can be either an integer or a vector of integers with components 
% increasing in size. In the first case, PCG is applied to the problem 
% with an mxm mesh. In the second case, PCG is applied to problems with 
% m(j)xm(j) meshes, j = 1, ..., length(m). 
%
% NOTE: If m is an integer vector, the components must be increasing 
% in size. 
%
% The inputs are follows:
% 
% Required:
%   m         = integer scalar or integer vector determing the number(s) of 
%               interior mesh points per side. 
%   Pre_Name  = string specifying name of selected preconditioner; 
%               acceptable values are 'none', 'Jacobi', 'block Jacobi', 
%               'SSOR', 'incomplete Cholesky', and -- for the Poisson 
%		problem ONLY -- 'ADI'. 
%   Problem   = string specifying name of the problem; acceptable values 
%               are 'Poisson' and 'biharmonic'. 
%   fun       = function handle (@functionname) for the right-hand 
%	        side f. 
%
% Optional:
%   tol       = stopping tolerance on the residual norm 
%	        (default 1.e-4*norm(initial right-hand side)). 
%   itmax     = maximum allowable number of iterations 
%	        (default 5000). 
%
% NOTES: (1) The demo uses the MATLAB pcg routine. 
%        (2) In block-Jacobi preconditioning, the block size is 
%	     currently determined automatically. This can easily 
%	     be changed to a user-specified value. 
%	 (3) SSOR preconditioning is currently used with omega = 1. 
%	     This can easily be changed to a user-specified value. 
%        (4) Incomplete Cholesky preconditioning currently allows 
%	     no fill-in. This can easily be changed. 
%        (5) In ADI preconditioning, the choice of rho is highly 
%	     unsophisticated. 
%
% For each value of m, the demo prints out the time (in seconds), the 
% number of iterations, the final relative residual norm, and an output 
% flag indicating success or failure (see the MATLAB on-line help for 
% flag details). The demo also plots the log of the relative residual 
% norms (the residual norms divided by the initial residual norm) vs. 
% the number of iterations for each mesh specified by m. 
%
% The initial guess for each PCG run is x = 0. 
%
% M-files required: 
%   adi_pre.m 
%   biharmonic.m
%   fun.m (called by loadrhs.m) 
%   laplacian.m 
%   loadrhs.m 
% 
if nargin < 4
   error('There must be at least four arguments (m,Pre_Name,Problem,fun).');
elseif nargin >  6
   error('There must be no more than six arguments.'); 
end
if nargin < 6 
   itmax = 5000; 
end 
global H_GL V_GL; % (for use in ADI preconditioning)
if ~(strcmp(Pre_Name,'none')|strcmp(Pre_Name,'Jacobi')|...
     strcmp(Pre_Name,'block Jacobi')|strcmp(Pre_Name,'SSOR')|...
     strcmp(Pre_Name,'incomplete Cholesky')|strcmp(Pre_Name,'ADI')) 
  error('Preconditioner name is not an acceptable string; see the comments.');
end
if ~(strcmp(Problem,'Poisson')|strcmp(Problem,'poisson')|...
     strcmp(Problem,'biharmonic'))
  error('Problem name is not an acceptable string; see the comments.');
end
if strcmp(Problem,'biharmonic')&strcmp(Pre_Name,'ADI') 
  error ('Cannot use ADI preconditioning on the biharmonic problem.'); 
end
set(figure(1),'Units','normal'); clf;
figure(1); hold on  
for j = length(m):-1:1 
    x = zeros(m(j)^2,1);
    if strcmp(Problem,'Poisson')|strcmp(Problem,'poisson')
      [A,H_GL,V_GL] = laplacian(m(j));
      A = -A; H_GL = -H_GL; V_GL = - V_GL;
    end
    if strcmp(Problem,'biharmonic')
      A = biharmonic_op(m(j));
    end
    b = loadrhs(fun,m(j)); 
    if nargin < 5
       tol = 1.e-4*((m(j)\1)*norm(b));
    end 
    tic; 
    if strcmp(Pre_Name,'none') 
        [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax); 
    end
    if strcmp(Pre_Name,'Jacobi') 
       M = diag(diag(A)); 
       [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,M);
    end
    if strcmp(Pre_Name,'block Jacobi') 
      n = m(j)^2; 
      k = fix(n/m(j)); 
      if k*m(j) ~= n
        error('The order of A must be divisible by the block size.');
      end
      if issparse(A)
        M = sparse(n,n); 
      else
        M = zeros(n,n);
      end
      for i = 1:k 
        M((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)) = ... 
                          A((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)); 
      end
      R = chol(M);
      [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,R',R);
    end
    if strcmp(Pre_Name,'SSOR')
      L = tril(A);
      D = diag(diag(A)); 
      R = D\(L');
      [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,L,R);
    end
    if strcmp(Pre_Name,'incomplete Cholesky')
      R = cholinc(A,'0');
      [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,R',R);
    end
    if strcmp(Pre_Name,'ADI')
      [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,@adi_pre);
    end
    time = toc;
    itno = size(resvec,1)-1; 
    logrelresvec = log10(resvec(1,1)\resvec);
    plot([0:1:itno]',logrelresvec); 
    text(itno,logrelresvec(size(logrelresvec,1),1)+.1,['m=' num2str(m(j))]); 
    drawnow; resvec = 0; 
    fprintf('\n For m = %d:  time = %g sec., its. = %d, relres = %g, flag = %d', ... 
                m(j), time, iter, relres, flag); 
    if j == length(m)
       y = x; 
    end
end 
fprintf('\n\n');
title(['PCG, preconditioner = ',Pre_Name,', ',Problem, ' problem']); 
xlabel('iteration number'); 
ylabel('log_{10} relative residual norm'); 
set(figure(1),'Units','pixels');
hold off
