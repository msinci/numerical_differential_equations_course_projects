function x = newtonbt(F,Fpr,x_0,tol_F,tol_x,itmax)
%
% Implements Newton's method with backtracking to solve a system 
% F(x) = 0. 
%
% The names of M-files for evaluating F(x) and F'(x) for a 
% given x are passed in as function handles (see below). 
%
% The user inputs one or more column vectors of starting values. 
% For each, the Newton iterations proceed until either the 
% residual norm is < tol_F, the step norm is < tol_x, or the 
% number of iterations is >= itmax. 
%
% The inputs are
% Required:
%   F        = function handle (@F) specifying M-file for 
%	       evaluating F(x)
%   Fpr      = function handle (@Fpr) specifying M-file for 
%	       evaluating F'(x)
%   x_0      = nx1 starting vector or nxm matrix, the columns of which 
%              are starting vectors. 
% Optional:
%   tol_F    = stopping tolerance: norm F(x) < tol_F => stop
%   tol_x    = stopping tolerance: step norm < tol_x => stop
%   itmax    = stopping tolerance: itno > itmax => stop
%
% The output is 
%   x        = final approximate solution 
%

% Check inputs and set the itmax stopping tolerance. 

if nargin > 6, error('Too many arguments.'); end
if nargin < 3, error('Not enough arguments.'); end 
if nargin < 6, itmax = 100; end   

t = 1.e-4;
thetamax = 1/2;
thetamin = 1/10;
nbtmax = 10;

% Loop over the number of starting values.

for run = 1:size(x_0,2)
  x = x_0(:,run);
  F_val = feval(F,x);
  F_norm = norm(F_val);
  itno = 0; 

% Set the remaining stopping tolerances. 

  if nargin < 5, 
     tol_x = (1+norm(x))\sqrt(eps); 
     if nargin < 4, 
       tol_F = (1+F_norm)\sqrt(eps); 
     end 
  end 

%>>> Print the output header and initial data. 
fprintf('\nIt.No. \t ||F(x)|| \t NB \t Step Red. Factor \n');
fprintf('%d \t %e\n', itno,F_norm);
%<<<

% Top of the Newton iteration loop. 

  while F_norm >= tol_F & itno < itmax

% Evaluate Fpr, the initial trial step, and F(x+step).  

    Fpr_val = feval(Fpr,x);        % For analytic J-evaluation. 
%  Fpr_val = feval(Fpr,x,F_val,F);  % For finite-difference J-evaluation. 
    if ~issparse(Fpr_val) 
      if cond(Fpr_val) > 1/(10*eps) 
          error('Ill-conditioned Jacobian.');
      end 
    end
    step = -Fpr_val\F_val;
    x_plus = x + step;
    F_val_plus = feval(F,x_plus);

% Set up for testing and possible backtracking. 

    lambda = 1; 
    nbt = 0;
    rho = norm(F_val_plus)/F_norm;

% Test the step and backtrack as necessary. 

    while rho > 1-t*lambda 
      nbt = nbt + 1;
      if nbt > nbtmax, error('Maximum number of backtracks reached.'); end 
      delta = rho^2 - 1 + 2*lambda;
      if delta <= 0 
        theta = thetamax;
      else
        theta = lambda/delta;
	theta = min(theta,thetamax);
	theta = max(theta,thetamin); 
      end
      step = theta*step;
      lambda = theta*lambda;
      x_plus = x + step;
      F_val_plus = feval(F,x_plus);
      rho = norm(F_val_plus)/F_norm;
    end

% Successful step. Update for the next Newton iteration. 

    x = x_plus;
    F_val = feval(F,x);
    F_norm = norm(F_val);
    itno = itno + 1;

%>>> Print iteration data. 
fprintf('%d \t %e \t %d \t %g \n', itno,F_norm,nbt,lambda);
%<<<

    step_norm = norm(step);
    if step_norm < tol_x, break, end

% End of Newton iteration loop. 

  end

  if F_norm <= tol_F,
    fprintf('Successful termination with ||F(x)|| < tol_F = %e.\n',tol_F); 
  elseif step_norm < tol_x,
    fprintf('Termination with ||step|| < tol_x = %e.\n',tol_x); 
  else 
   fprintf('Termination with iteration no. = itmax = %d.\n',itmax); 
  end

% End of loop over the starting values. 

end
