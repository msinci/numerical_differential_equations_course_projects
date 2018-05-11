% Demonstrates MATLAB's ode45 on the HW2 IVP 
%
% y' = y + exp(t)*cos(t), y(0) = 0.
%
% The solution is y(t) = exp(t)*sin(t).

% First, use ode45 defaults (i.e., no options). 

y_0 = 0;
t_0 = 0; 
t_f = pi;

fun = @(t,y)y + exp(t)*cos(t);

% Call ode45.
[tvals,yvals] = ode45(fun,[t_0 t_f],y_0);

maxerr = max(abs(yvals - exp(tvals).*sin(tvals)));
nsteps = size(tvals,1);

fprintf('\n nsteps = %d, no. f-evals = %d, maxerr = %g \n',...
    nsteps, 6*nsteps, maxerr);

figure(1); plot(tvals,yvals);

%% Now illustrate setting tolerances using ode45 options 

y_0 = 0;
t_0 = 0; 
t_f = pi;

fun = @(t,y)y + exp(t)*cos(t);

% Set the error tolerances to 1e-8 in structure "myopts". 
myopts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% Call ode45.
[tvals,yvals] = ode45(fun,[t_0 t_f],y_0,myopts);

maxerr = max(abs(yvals - exp(tvals).*sin(tvals)));
nsteps = size(tvals,1);

fprintf('\n nsteps = %d, no. f-evals = %d, maxerr = %g \n',...
    nsteps, 6*nsteps, maxerr);

figure(1); plot(tvals,yvals);

%% Now illustrate event location using ode45 options.  
% This sets the options structure to terminate the integration
% when the solution reaches a maximum value, signaled by an  
% "event" function that detects when the derivative goes from 
% positive to negative. 
y_0 = 0;
t_0 = 0; 
t_f = pi;

fun = @(t,y)y + exp(t)*cos(t);

% Set the options structure to call the event function, defined  
% in a separate m-file because it must return multiple values. 

myopts = odeset('Events',@hw2funevent);

% Call ode45 with the options structure as the last argument.
[tvals,yvals,tEvent,yEvent] = ode45(fun,[t_0 t_f],y_0,myopts);

maxerr = max(abs(yvals - exp(tvals).*sin(tvals)));
nsteps = size(tvals,1);

fprintf('\n nsteps = %d, no. f-evals = %d, maxerr = %g \n',...
    nsteps, 6*nsteps, maxerr);
yf = yvals(nsteps); tf = tvals(nsteps);
fprintf('\n Max value of y = %g; occurs at t = %g \n',yf,tf);

figure(1); plot(tvals,yvals);