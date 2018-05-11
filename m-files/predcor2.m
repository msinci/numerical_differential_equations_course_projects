function [T,Y,nsteps,nfe,nje] = predcor2(odefun,t_int, y_0, itmeth, odejac) 
%
% Implements the second-order implicit Adams-Moulton 
% (trapezoidal) method with variable stepsizes using the 
% second-order explicit Adams-Bashforth method as a predictor. 
% Stepsizes are determined to maintain a bound on local error per 
% unit step. Local error is estimated using the difference between 
% the "predicted" Adams-Bashforth point and the Adams-Moulton 
% (trapezoidal) point at each step. The forward Euler and backward 
% Euler methods are used for startup, with local error estimated 
% using the difference in the forward- and backward-Euler points. 
%  
% The inputs are as follows: 
%
% Required:
%   odefun     = function handle (@odefun) for evaluating the 
%	         right-hand side of the ODE. This should 
%		 accept (t,y) as arguments and return the value 
%		 of the function. 
%   t_int      = two-vector [t_0, t_f] of initial and final times.
%   y_0	       = initial value of the solution (column vector). 
%
% Optional: 
%   itmeth     = string specifying the type of corrector iteration. 
%		 Must be either 'fixed-point' (default) or 'Newton'. 
%		 If itmeth = 'Newton', then the odejac argument must 
%		 be provided. 
%   odejac     = function handle (@odejac) for evaluating the 
%		 Jacobian with respect to y of the right-hand 
%		 side of the ODE. This should accept (t,y) as 
%		 arguments and return the value of the Jacobian. 
%		 Required when itmeth = 'Newton'. 
%
% The outputs are as follows: 
%
%   T          = (nsteps + 1)-vector of equally spaced t-values. 
%   Y	       = (nsteps + 1)xn matrix, where n is the length
%		 of y_0. The ith row of Y is the transpose of 
%		 the approximate solution at the ith value of T. 
%   nsteps     = number of steps
%   nfe	       = number of function evaluations
%   nje	       = number of Jacobian evaluations 
%

% Check inputs and set defaults.
if ~isa(odefun,'function_handle'),
   error('The odefun argument must be a function handle.');
end 
if size(t_int) ~= [1 2], 
  error('t_int must be a 1x2 vector [t_0,t_f] with t_0 < t_f.'); 
end 
if t_int(1,1) >= t_int(1,2), 
  error('t_int must be a 1x2 vector [t_0,t_f] with t_0 < t_f.'); 
end
if size(y_0,1) < size(y_0,2), error('y_0 must be a column vector.'); end
if nargin < 3, error('Not enough arguments.'),end
if nargin < 4, itmeth = 'fixed-point'; end 
if ~strcmp(itmeth,'fixed-point') & ~strcmp(itmeth,'Newton'), 
  error('The itmeth flag must be either ''fixed-point'' or ''Newton''.');
end 
if strcmp(itmeth,'Newton'), 
  if nargin < 5, 
    error('The odejac argument must be provided when itmeth = ''Newton''.');
  end
  if ~isa(odejac,'function_handle'),
    error('The odefun argument must be a function handle.');
  end 
end

% Initialize tolerances and counters. 
abstol = 1.e-6;                   % Absolute local error tolerance.
reltol = 1.e-5;			  % Relative local error tolerance. 
itreltol  = 1.e-4;		  % Relative corrector-iteration tolerance. 
if strcmp(itmeth,'fixed-point'),  % Max number of corrector iterations ... 
   itmax = 5;			  % ... for fixed-point iteration ...
else
   itmax = 3;			  % ... and for Newton iteration. 
end

maxsteps = 10000;                 % Maximum number of time steps.
nsteps = 0;                       % Time step counter.
nfe = 0;			  % Function-evaluation counter. 
nje = 0;			  % Jacobian-evaluation counter. 

% Initialize problem-specific quantities. 
t_0 = t_int(1,1);                 % Initial t-value. 
t_f = t_int(1,2);		  % Final t-value. 
h = 1000\(t_f - t_0);		  % Initial stepsize. 
hmin = 10000\h;			  % Minimum allowable stepsize. 

n = size(y_0,1); 
T = t_0;
Y = y_0'; 

tcur = t_0; 
ycur = y_0; 
fcur = feval(odefun,tcur,ycur); 
nfe = nfe + 1;
Y = ycur'; 
T = tcur;
hpast = h;

% Top of the integration loop. 
while tcur < t_f, 
  if nsteps > maxsteps, error('Maximum number of steps reached.'); end 
% Check the stepsize. Adjust if necessary in order not to 
% overshoot t_f.
  h = min(h,t_f - tcur);
  errtol = abstol + reltol*norm(ycur);
  stepflag = 0;
  if nsteps < 2 & strcmp(itmeth,'Newton'), hjac = h; jflag = 1; end
  while stepflag == 0, 
    tnext = tcur + h;
%   Reset the past f-value if h has changed (after the first step). 
    if h ~= hpast & nsteps > 0, 
      fpast = fcur + (hpast\h)*(fpast - fcur);
      hpast = h; 
    end
%   Compute the predictor step using ... 
    if nsteps == 0,
%   ... forward Euler ... 
      ynextpred = ycur + h*fcur; 
    else
%   ... or 2nd order Adams-Bashforth.
      ynextpred = ycur + (2\h)*(3*fcur - fpast); 
    end
    ynext = ynextpred;
    fnext = feval(odefun,tnext,ynext); 
    nfe = nfe + 1;
    if nsteps == 0, 
      residnorm = norm(ynext - (ycur + h*fnext)); 
    else 
      residnorm = norm(ynext - (ycur + (2\h)*(fcur + fnext))); 
    end
    ittol = itreltol*residnorm; 
    itno = 0; 
%   Implement the corrector iterations ... 
    while residnorm > ittol & itno <= itmax, 
      itno = itno + 1; 
      if itno > itmax, break; end
      if nsteps == 0, 
%     ... for backward Euler ... 
        if strcmp(itmeth,'fixed-point'),
          ynext = ycur + h*fnext; 
        else
          if jflag == 1 | hjac\abs(h-hjac) > .3, 
            Jval = eye(n,n) - h*feval(odejac,tnext,ynext);
	    nje = nje + 1;
	    jflag = 0;
	    hjac = h;
	    [L,U] = lu(Jval);
	  end 
          ynext = ynext - U\(L\(ynext - (ycur + h*fnext))); 
        end
      else
%     ... or the trapezoidal method. 
        if strcmp(itmeth,'fixed-point'),
          ynext = ycur + (2\h)*(fcur + fnext); 
        else
          if jflag == 1 | hjac\abs(h-hjac) > .3, 
	    Jval = eye(n,n) - (2\h)*feval(odejac,tnext,ynext);
	    nje = nje + 1;
	    jflag = 0;
	    hjac = h;
	    [L,U] = lu(Jval);
	  end 
          ynext = ynext - U\(L\(ynext - (ycur + (2\h)*(fcur + fnext)))); 
        end
      end 
      fnext = feval(odefun,tnext,ynext); 
      nfe = nfe + 1;
      residnormold = residnorm; 
      if nsteps == 0, 
        residnorm = norm(ynext - (ycur + h*fnext));
      else
        residnorm = norm(ynext - (ycur + (2\h)*(fcur + fnext))); 
     end 
      if residnorm > residnormold, break; end
    end
%   Estimate the local error.
    if nsteps == 0, 
      errest = 2\norm(ynext - ynextpred); 
      if errest > 0,
        hplus = ((errest\.9*errtol)^(2\1))*h;
      else 
        hplus = 4*h;
      end 
     else
      errest = 6\norm(ynext - ynextpred); 
      if errest > 0,
        hplus = ((errest\.9*errtol)^(3\1))*h;
      else 
        hplus = 4*h;
      end 
    end 
%   If the error is too large or the corrector iterations failed, 
%   try again to determine a successful step. 
    if errest > errtol | itno > itmax | residnorm > residnormold,  
%     If using Newton's method and the Jacobian is not fresh, 
%     refresh it and try the step again. 
      if strcmp(itmeth,'Newton')&(itno>itmax|residnorm>residnormold)&jflag==-1,
        jflag = 1;
      else
%     Otherwise, try a reduced step. 
        h = min(hplus,2\h);
	jflag = -1;
        if h < hmin, error('Minimum stepsize reached.'); end
      end
    else
%   If the error is small enough, accept the step and adjust h 
%   for the next step. 
      stepflag = 1;
      jflag = -1;
      h = min(hplus,4*h); 
      nsteps = nsteps + 1;
    end
  end
% Update the Y-array and other variables. 
  Y = [Y;ynext']; 
  if nsteps == 1, hpast = h; end
  T = [T;tnext];  
  tcur = tnext;
  ycur = ynext; 
  fpast = fcur;
  fcur = fnext;
end
% End of the integration loop. 
