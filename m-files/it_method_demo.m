function it_method_demo(m,method,fun,tol,itmax)
%
% Demonstrates the performance of Jacobi iteration, the Gauss-Seidel 
% method, or SOR on the Poisson problem 
%
%             (Laplacian)u = f   in D
%                        u = 0   on the boundary of D
%
% where D is the unit square in R^2 discretized on an evenly 
% spaced mxm mesh of interior points. 
%
% Specifically, Jacobi, Gauss-Seidel, or SOR (with the optimal 
% acceleration parameter) is applied to the problem for a given 
% f and square meshes with size determined by m. 
%
% m can be either an integer or a vector of integers with components 
% increasing in size. In the first case, Jacobi, Gauss-Seidel, or SOR
% is applied to the problem with an mXm mesh. In the second case, 
% Jacobi, Gauss-Seidel, or SOR is applied to problems with m(j)Xm(j) 
% meshes, j = 1, ..., length(m). 
%
% NOTE: If m is an integer vector, the components must be increasing 
% in size. 
%
% The inputs are as follows:
% 
% Required:
%   m       = scalar integer or integer vector determining the 
%             number(s) of interior mesh points per side. 
%   method  = 'Jacobi' or 'Gauss-Seidel' or 'SOR' (string) 
%   fun     = function handle (@functionname) for the right-hand 
%	      side f. 
% Optional:
%   tol     = stopping tolerance on the residual norm 
%	      (default 1.e-4*norm(initial right-hand side)). 
%   itmax   = maximum allowable number of iterations 
%	      (default 5000). 
%
% The output is a plot of the log of the relative residual norms 
% (the residual norms divided by the initial residual norm) vs. 
% the number of iterations for each mesh specified by m. 
%
% The initial guess for each Jacobi or Gauss-Seidel run is x = 0. 
%
% M-files required: 
%   laplacian.m              (sets up the discretized Laplacian)
%   jacobi.m		     (performs Jacobi iteration)
%   gaussseidel.m	     (performs Gauss-Seidel iteration)
%   sor.m		     (performs SOR) 
%   loadrhs.m		     (loads the right-hand side f)
% 
%
set(figure(1),'Units','normal');clf;
figure(1); hold on; 
if nargin < 3
   error('There must be at least three arguments (m, method, fun).');
elseif nargin >  5
   error('There must be no more than five arguments.');
end
if nargin < 5 
   itmax = 5000; 
end 
for j = length(m):-1:1 
    x = zeros(m(j)^2,1);
    A = laplacian(m(j));
    b = loadrhs(fun,m(j));
    if nargin < 4
       tol = 1.e-4*norm(b);
    end 
    tic; 
    if (strcmp(method,'Jacobi')|strcmp(method,'jacobi'))
      [x,itno,resvec] = jacobi(x, A, b, tol, itmax);
    elseif (strcmp(method,'Gauss-Seidel')|strcmp(method,'gauss-seidel'))
      [x,itno,resvec] = gaussseidel(x, A, b, tol, itmax);
    elseif (strcmp(method,'SOR')|strcmp(method,'sor'))
      w = 2/(1 + sin(pi/(m(j)+1)));
      [x,itno,resvec] = sor(x, A, b, w, tol, itmax);
    else
      error('method must be either Gauss-Seidel, Jacobi, or SOR in single quotes');
    end
    time = toc;
    logrelresvec = log10(resvec(1,1)\resvec);
    plot(itno,logrelresvec);
    text(itno(size(itno,2))-1,logrelresvec(size(logrelresvec,2))-.1,['m=' num2str(m(j))]);
    drawnow 
    fprintf('\n For m = %d, run time = %g seconds', m(j), time);
    if j == length(m)
       y = x; 
    end
end 
fprintf('\n\n');
if (strcmp(method,'Gauss-Seidel')|strcmp(method,'gauss-seidel')) 
  title('Gauss-Seidel iteration, Poisson problem' ); 
elseif (strcmp(method,'SOR')|strcmp(method,'sor'))
  title('SOR, Poisson problem');
else
  title('Jacobi iteration, Poisson problem'); 
end
xlabel('iteration number'); 
ylabel('log_{10} relative residual norm'); 
set(figure(1),'Units','pixels');
hold off

