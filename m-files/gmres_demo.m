function gmres_demo(m,c,d,e,restart,Pre_Name,fun,tol,itmax)
%
% Demonstrates the performance of the generalized minimal residual 
% method (GMRES) on the problem 
%
%           -(Laplacian)u + cu + du_x + eu_y = f 
%
% with Dirichlet boundary conditions on the unit square in R^2 
% discretized on an evenly spaced mxm mesh of interior points. 
%
% m can be either an integer or a vector of integers with components 
% increasing in size. In the first case, GMRES is applied to the problem 
% with an mxm mesh. In the second case, GMRES is applied to problems with 
% m(j)xm(j) meshes, j = 1, ..., length(m). 
%
% NOTE: If m is an integer vector, the components must be increasing 
% in size. 
%
% The inputs are as follows:
% 
% Required:
%   m         = integer scalar or integer vector determing the number(s) of 
%               interior mesh points per side.
%   c         = coefficient of zero-order term. 
%   d         = coefficient of x-derivative term. 
%   e         = coefficient of y-derivative term. 
%   restart   = number of GMRES steps before restarting 
%   Pre_Name  = string specifying name of selected preconditioner; 
%               acceptable values are 'none', 'Jacobi', 'block Jacobi', 
%               'ILU', or 'ADI'. 
%   fun       = function handle (@functionname) for the right-hand 
%	        side f. 
%
% Optional:
%   tol       = stopping tolerance on the residual norm 
%		(default 1.e-4*norm(initial right-hand side)). 
%   itmax     = maximum allowable number of iterations 
%		(default 5000). 
%
% NOTES: (1) The demo uses the MATLAB gmres routine. 
%        (2) In block-Jacobi preconditioning, the block size is 
%	     currently determined automatically. This can easily 
%	     be changed to a user-specified value. 
%        (3) ILU preconditioning currently allows no fill-in. 
%	     This can easily be changed. 
%        (4) In ADI preconditioning, the choice of rho is highly 
%	     unsophisticated. 
%
% For each value of m, the demo prints out the time (in seconds), the 
% number of iterations, the final relative residual norm, and an output 
% flag indicating success or failure (see the MATLAB on-line help for 
% flag details). The demo also plots the log of the relative residual 
% norms (the residual norms divided by the initial residual norm) vs. 
% the number of iterations for each mesh specified by m. 
%
% The initial guess for each GMRES run is x = 0. 
%
% M-files required: 
%   adi_pre.m 
%   condif.m
%   fun.m (called by loadrhs.m) 
%   loadrhs.m 
% 
if nargin < 7
   error('There must be at least seven arguments (m,c,d,e,restart,Pre_Name,fun).');
elseif nargin >  9
   error('There must be no more than nine arguments.'); 
end
if nargin < 9 
   itmax = 5000; 
end 
maxit = floor(restart\itmax); 
global H_GL V_GL; % (for use in ADI preconditioning)
if ~(strcmp(Pre_Name,'none')|strcmp(Pre_Name,'Jacobi')|...
     strcmp(Pre_Name,'block Jacobi')|...
     strcmp(Pre_Name,'ILU')|strcmp(Pre_Name,'ADI')) 
  error('Preconditioner name is not an acceptable string. See comments.');
end
fprintf('\n\n GMRES(%d), preconditioner = %s',restart,Pre_Name);
fprintf('\n Poisson problem, c = %g, d = %g, e = %g', c, d, e); 
set(figure(1),'Units','normal'); clf;
figure(1); hold on; 
for j = length(m):-1:1 
    x = zeros(m(j)^2,1);
    [A,H_GL,V_GL] = condif(c,d,e,m(j));
    b = loadrhs(fun,m(j)); 
    tic; 
    if strcmp(Pre_Name,'none') 
       if nargin < 8
          tol = 1.e-4*((m(j)\1)*norm(b));
       end 
       [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit); 
    end
    if strcmp(Pre_Name,'Jacobi') 
       M = diag(diag(A)); 
       if nargin < 8
          tol = 1.e-4*((m(j)\1)*norm(M\b));
       end 
       [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,M);
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
      [M1,M2] = lu(M); 
       if nargin < 8
          tol = 1.e-4*((m(j)\1)*norm(M2\(M1\b)));
       end 
      [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,M1,M2);
    end
    if strcmp(Pre_Name,'ILU')
       [M1,M2] = luinc(A,'0'); 
       if nargin < 8
          tol = 1.e-4*((m(j)\1)*norm(M2\(M1\b)));
       end 
       [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,M1,M2);
    end
    if strcmp(Pre_Name,'ADI')
       if nargin < 8
          tol = 1.e-4*((m(j)\1)*norm(adi_pre(b)));
       end 
       [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,@adi_pre);
    end
    time = toc;
    itno = size(resvec,1)-1; 
    logrelresvec = log10(resvec(1,1)\resvec);
    plot([0:1:itno]',logrelresvec); 
    text(itno,logrelresvec(size(logrelresvec,1),1)+.1,['m=' num2str(m(j))]); 
    drawnow; resvec = 0; 
    fprintf('\n For m = %d:  time = %g sec., its. = %d, relres = %g, flag = %d', ... 
                m(j), time, restart*iter(1)+iter(2), relres, flag); 
    if j == length(m)
       y = x; 
    end
end 
fprintf('\n\n'); 
title(['GMRES(',num2str(restart),'), preconditioner = ',Pre_Name, ', Poisson problem']);
xlabel('iteration number'); 
ylabel('log_{10} relative residual norm'); 
set(figure(1),'Units','pixels');
hold off


