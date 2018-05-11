function pcg_demo(m,PreName,Problem,RHSfun,tol,itmax)
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
%   PreName  = string specifying name of selected preconditioner; 
%               acceptable values are 'none', 'Jacobi', 'block Jacobi', 
%               'SSOR', 'incomplete Cholesky', and -- for the Poisson 
%		problem ONLY -- 'ADI'. 
%   Problem   = string specifying name of the problem; acceptable values 
%               are 'Poisson' and 'biharmonic'. 
%   RHSfun    = function handle (@RHSfunctionname) for the right-hand 
%	            side f. 
%
% Optional:
%   tol       = stopping tolerance on the relative residual norm, i.e.,  
%             the residual norm divided by the initial residual norm. 
%	          (default 1.e-6). 
%   itmax     = maximum allowable number of iterations (default 5000).
%
% NOTES: (1) The demo uses the MATLAB pcg routine. 
%        (2) In block-Jacobi preconditioning, the block size is 
%	     currently determined automatically. This can easily 
%	     be changed to a user-specified value. 
%	     (3) SSOR preconditioning is currently used with omega = 1. 
%	     This can easily be changed to a user-specified value. 
%        (4) Incomplete Cholesky preconditioning currently allows 
%	     no fill-in. This can easily be changed. 
%        (5) In ADI preconditioning, the choice of rho is naive; other 
%        choices may perform better. 
%
% For each value of m, the demo prints out the time (in seconds), the 
% number of iterations, the final relative residual norm (the residual
% norm divided by the initial residual norm), and an output  
% flag indicating success or failure (see the MATLAB on-line help for 
% flag details). The demo also plots the log of the relative residual 
% norms  vs. 
% the number of iterations for each mesh specified by m. 
%
% The initial guess for each PCG run is x = 0. 
%
% M-files required: 
%   adi_pre.m     (applies ADI preconditioning)
%   biharmonic.m  (sets up the discretized biharmonic operator)
%   RHSfun.m      (user-supplied, called by loadrhs.m) 
%   laplacian.m   (sets up the discretized laplacian)
%   loadrhs.m     (loads the right-hand side)
% 

global H_GL V_GL; % (for use in ADI preconditioning)

if nargin < 4
   error('There must be at least four arguments (m,PreName,Problem,RHSfun).');
elseif nargin >  6
   error('There must be no more than six arguments.'); 
end
if nargin < 6 
   itmax = 5000; 
end 
if nargin < 5
   tol = 1.e-6;
end 

if ~(strcmp(PreName,'none')|strcmp(PreName,'Jacobi')|...
     strcmp(PreName,'block Jacobi')|strcmp(PreName,'SSOR')|...
     strcmp(PreName,'incomplete Cholesky')|strcmp(PreName,'ADI')) 
  error('Preconditioner name is not an acceptable string; see the comments.');
end
if ~(strcmp(Problem,'Poisson')|strcmp(Problem,'poisson')|...
     strcmp(Problem,'biharmonic'))
  error('Problem name is not an acceptable string; see the comments.');
end
if strcmp(Problem,'biharmonic')&strcmp(PreName,'ADI') 
  error ('Cannot use ADI preconditioning on the biharmonic problem yet.'); 
end
set(figure(gcf),'Units','normal'); clf;
figure(gcf); hold on  

fprintf(['\n PCG, Preconditioner = ' PreName]);

for j = length(m):-1:1 
    if strcmp(Problem,'Poisson')|strcmp(Problem,'poisson')
      [A,H_GL,V_GL] = laplacian(m(j));
      A = -A; H_GL = -H_GL; V_GL = - V_GL;
    end
    if strcmp(Problem,'biharmonic')
      A = biharmonic_op(m(j));
    end
    b = loadrhs(RHSfun,m(j)); 
    tic; 
    if strcmp(PreName,'none') 
        M = [];
    end
    if strcmp(PreName,'Jacobi') 
       M = diag(diag(A)); 
    end
    if strcmp(PreName,'block Jacobi') 
      n = m(j)^2; 
      k = fix(n/m(j)); 
      if k*m(j) ~= n
        error('The order of A must be divisible by the block size.');
      end
      if issparse(A)
        MM = sparse(n,n); 
      else
        MM = zeros(n,n);
      end
      for i = 1:k 
        MM((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)) = ... 
                          A((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)); 
      end
      R = chol(MM);
      M = @(v)R\(R'\v);
    end
    if strcmp(PreName,'SSOR')
      omega = (1 + sin((m(j)+1)\pi))\2;
      D = omega\diag(diag(A));
      L = D + tril(A,-1);
      M = @(v)(2-omega)*(L'\(D*(L\v)));
    end
    if strcmp(PreName,'incomplete Cholesky')
%        opts.type = 'nofill';
      R = ichol(A,struct('type','nofill'));
      M = @(v)R\(R'\v);
    end
    if strcmp(PreName,'ADI')
      M = @adi_pre;
    end
    [x,flag,relres,iter,resvec] = pcg(A,b,tol,itmax,M);
    time = toc;
    itno = size(resvec,1)-1; 
    logrelresvec = log10(resvec(1,1)\resvec);
    figure(gcf);plot((0:1:itno)',logrelresvec); 
    text(itno,logrelresvec(size(logrelresvec,1),1)+.1,['m=' num2str(m(j))]); 
    drawnow; resvec = 0; 
    fprintf('\n For m = %d:  time = %g sec., its. = %d, relres = %g, flag = %d', ... 
                m(j), time, iter, relres, flag); 
    if j == length(m)
       y = x; 
    end
end 
fprintf('\n\n');
title(['PCG, preconditioner = ',PreName,', ',Problem, ' problem']); 
xlabel('iteration number'); 
ylabel('log_{10} relative residual norm'); 
set(figure(gcf),'Units','pixels');
hold off
