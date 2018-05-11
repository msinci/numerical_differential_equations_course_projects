function x = gmres_demo(m,coeffs,RHSfun,precond,tol,restart,maxrestart)
%
% Demonstrates the performance of MATLAB's GMRES with optional left
% preconditioning with ILU (incomplete LU factorization) on the
% problem  
%
%             -(Laplacian) u + cu + du_x + eu_y = f  
%
% with Dirichlet boundary conditions on the unit square in R^2 
% discretized on an evenly spaced mxm mesh of interior points in
% the unit square [0,1]x[0,1]. 
%
% m can be either an integer or a vector of integers with components 
% increasing in size. In the first case, GMRES is applied to the problem 
% with an mxm mesh. In the second case, GMRES is applied to
% problems with  m(j)xm(j) meshes,  j = 1, ..., length(m). 
%
% NOTE: If m is an integer vector, make sure the components are 
% increasing in size. 
%
% Currently, the preconditioning options are 'none' (no
% preconditioning) and 'ILU' (incomplete LU using MATLAB's default settings).
%
% The inputs are as follows: 
% Required:
%   m         = integer scalar or integer vector determining the 
%               number(s) of interior mesh points per side. 
%   coeffs    = row vector of coefficients [c d e].
%   RHSfun    = function handle (@functionname) for evaluating the 
%		right-hand side function f at given (x,y). The routine 
%		should accept scalar x and y as arguments and return f(x,y). 
%   precond   = string specifying the preconditioner; acceptable 
%               values are 'none', 'Jacobi', 'block Jacobi', 'ADI', and 'ILU'. 
% Optional:
%   tol       =  stopping tolerance on the residual norm (default 1.e-6). 
%   restart  = GMRES restart value (default 40).
%   maxrestart  = maximum allowable number of restarts (default 50). 
%
% NOTES: Jacobi preconditioning uses the matrix diagonal; block
% Jacobi uses mxm blocks on the matrix diagonal; ILU uses the
% MATLAB defaults (choosing various options in the MATLAB
% implementat would probably improve performance); ADI uses a naive
% parameter choice.  
%
% The output is a graph of iteration numbers versus log residual norms 
% for each mesh specified by m. The program also prints out the run time 
% (in seconds) for each value of m. 
%
% The initial guess for each GMRES run is x = 0. The iterations terminate 
% when the residual norm is <= tol*(initial residual norm); the
% default is tol = 1.e-6. For safety, GMRES is allowed to iterate
% for no more than restart*maxrestart iterations; with the default
% values, this is 40*50 = 2000. 
% 
% M-files required: 
%   adi_pre.m 
%   condif.m      (sets up the discretized convection-diffusion operator) 
%   RHSfun.m    (user-supplied, called by loadrhs.m) 
%   loadrhs.m    (loads the right-hand side)
%
% Optional M-file:
%   seesol.m     (shows the solution) 
% 

global H_GL V_GL; % (for use in ADI preconditioning)

if nargin < 4, 
  error('gmres_demo: The number of arguments must be => 4; see comments.'); 
elseif nargin > 7, 
  error('gmres_demo: The number of arguments must be <= 7; see comments.');
end

if ~(strcmp(precond,'none')|strcmp(precond,'Jacobi')|...
     strcmp(precond,'block Jacobi')|strcmp(precond,'ILU')|...
     strcmp(precond,'ADI'))
  error('gmres_demo: Preconditioner name is not an acceptable string; see comments.');
end

if nargin < 7
  maxrestart = 50;
elseif isempty(maxrestart)
  maxrestart = 50;
end

if nargin < 6, 
  restart = 40; 
elseif isempty(restart), 
  restart = 40; 
end 

if nargin < 5, 
  tol = 1.e-6; 
elseif isempty(tol), 
  tol = 1.e-6; 
end 

figure(gcf); clf; hold on  

fprintf(['\n GMRES(%d), Preconditioner = ' precond ], restart);

for j = length(m):-1:1 
  [A,H_GL,V_GL] =  condif(m(j),coeffs); 
  b = loadrhs(RHSfun,m(j)); 
  tic; 
    if strcmp(precond,'none')
      M = [];
    end
    if strcmp(precond,'Jacobi') 
       M = diag(diag(A)); 
    end
    if strcmp(precond,'block Jacobi') 
      n = m(j)^2; 
      if issparse(A)
        MM = sparse(n,n); 
      else
        MM = zeros(n,n);
      end
      for i = 1:m(j)
        MM((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)) = ... 
                          A((i-1)*m(j)+1:i*m(j),(i-1)*m(j)+1:i*m(j)); 
      end
      [L,U] = lu(MM); 
      M = @(v)U\(L\v);
    end
    if strcmp(precond,'ILU')
      [L,U] = ilu(A);
      M = @(v)U\(L\v);
    end
     if strcmp(precond,'ADI')
       M = @adi_pre;
    end
    
    [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxrestart,M);

     time = toc;

    if flag > 0
      fprintf('\n Warning: In gmres_demo, gmres flag = %d \n', flag),
    end

    lrrnrm = log10(resvec(1,1)\resvec);
    itno = 1:length(lrrnrm);

    figure(gcf); plot(itno-1,lrrnrm);
    text(itno(length(itno)),lrrnrm(length(lrrnrm))+.1,['m=' num2str(m(j))]);
    drawnow 

fprintf('\n For m = %d:  %d iterations; time = %g seconds', m(j), length(itno), time);
    if j == length(m)
       y = x; 
    end
end 
fprintf('\n\n');
title(['GMRES(',num2str(restart),'), preconditioner = ',precond,...
       ', Convection-Diffusion Problem: c = ',num2str(coeffs(1)),...
       ', d = ', num2str(coeffs(2)),', e = ', num2str(coeffs(3))]); 

xlabel('iteration number'); 
ylabel('log_{10} relative residual norm'); 
hold off
