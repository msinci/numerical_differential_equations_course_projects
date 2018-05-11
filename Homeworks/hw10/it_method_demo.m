function it_method_demo(m,method,RHSfun,tol,itmax)
%
% Demonstrates the performance of Jacobi iteration, the Gauss-Seidel 
% method, or SOR on the Poisson problem 
%
%            -(Laplacian)u = f   in D
%                        u = 0   on the boundary of D
%
% where D is the unit square in R^2 discretized on an evenly 
% spaced mxm mesh of interior points. 
%
% m can be either an integer or a vector of integers with components 
% increasing in size. In the first case, the method is applied to the 
% problem with an mXm mesh. In the second case, the method is applied 
% to problems with m(j)Xm(j) meshes, j = 1, ..., length(m). 
%
% NOTE: If m is an integer vector, the components must be increasing 
% in size. 
%
% The inputs are 
% Required:
%   m       = integer scalar or integer vector determining the 
%             number(s) of interior mesh points per side. 
%   method  = string specifying the method; acceptable values 
%             are 'Jacobi' or 'Gauss-Seidel' or 'SOR' (string). 
%   RHSfun   = function handle (@functionname) for evaluating the 
%              right-hand side function f at given (x,y). The routine 
%              should accept scalar x and y as arguments and return f(x,y). 
% Optional:
%   tol     = stopping tolerance on the relative residual norm 
%             (default 1.e-6). 
%   itmax   = maximum allowable number of iterations (default 3000). 
%
% The output is a graph of iteration numbers versus log residual norms 
% for each mesh specified by m. 
%
% The initial guess for each method run is x = 0. The iterations 
% terminate when the residual norm is less than or equal to tol; 
% the default is tol = 1.e-6*norm(b). For safety, the method is 
% allowed to iterate for no more than itmax iterations; the default 
% is itmax = 3000. 
% 
% M-files required: 
%   laplacian.m              (sets up the discretized Laplacian)
%   jacobi.m		     (performs Jacobi iteration)
%   gaussseidel.m	     (performs Gauss-Seidel iteration)
%   sor.m		             (performs SOR) 
%   RHSfun.m                (user-supplied, called by loadrhs.m) 
%   loadrhs.m		     (loads the right-hand side)
% 

if nargin < 3
  error('it_method_demo: The number of arguments must be => 3; see comments.'); 
elseif nargin >  5
  error('it_method_demo: The number of arguments must be <= 5; see comments.');
end

if nargin < 5, 
  itmax = 3000; 
elseif isempty(itmax), 
  itmax = 3000; 
end 

fig = gcf;

set(fig,'Units','normal');
clf

fprintf(['\n Method = ' method]);

for j = length(m):-1:1 
    x = zeros(m(j)^2,1);
    A = -laplacian(m(j));
    b = loadrhs(RHSfun,m(j));
    if nargin < 4, 
      tol = 1.e-6*norm(b); 
    elseif isempty(tol), 
      tol = 1.e-6*norm(b); 
    end 
    tic; 
    if (strcmp(method,'Gauss-Seidel')|strcmp(method,'gauss-seidel'))
      [x,itno,resvec] = gaussseidel(x, A, b, tol, itmax);
    elseif (strcmp(method,'Jacobi')|strcmp(method,'jacobi'))
      [x,itno,resvec] = jacobi(x, A, b, tol, itmax);
    elseif (strcmp(method,'SOR')|strcmp(method,'sor'))
      w = 2/(1 + sin(pi/(m(j)+1)));
      [x,itno,resvec] = sor(x, A, b, w, tol, itmax);
    else
      error('it_method_demo: The method must be either Gauss-Seidel, Jacobi, or SOR in single quotes;see comments');
    end
    time = toc;
     lrrnrm = log10(resvec(1,1)\resvec);
    figure(fig); 
   % plot(itno-1,lrrnrm - log10(m(j))); 
   plot(itno-1,lrrnrm); 
    text(itno(length(itno))-.1,lrrnrm(length(lrrnrm)),['m=' num2str(m(j))]);
    drawnow; figure(fig);
%    fprintf('\n For m = %d, run time = %g seconds', m(j), time);
    fprintf('\n For m = %d:  %d iterations; time = %g seconds', m(j), length(itno), time);
    hold on  
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
ylabel('log_{10} residual norm'); 
set(fig,'Units','pixels');
hold off

