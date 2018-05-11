%% Demonstrates forward Euler on the Gear IVP. 

y_0 = 0;
t_0 = 0; 
t_f = 1;
nsteps = 50000;

fun = @(t,y)-1000*(y-t^2) + 2*t;

[yvals,tvals] = forward_euler(fun,y_0,t_0,t_f,nsteps);

figure(1); plot(tvals,yvals);

maxerr = max(abs(yvals - tvals.^2));

fprintf('\n nsteps = %d, maxerr = %g \n',nsteps, maxerr);

%% Demonstrates backward Euler on the Gear IVP. 

y_0 = 0;
t_0 = 0; 
t_f = 1;
nsteps = 50;

h = nsteps\(t_f - t_0); 
tvals = [t_0:h:t_f]; 
yvals = zeros(size(y_0,1),nsteps + 1); 
yvals(:,1) = y_0; 
for n = 1:nsteps 
  yvals(:,n+1) = (1+1000*h)\ ...
      (yvals(:,n)+h*(1000*tvals(n+1)^2+2*tvals(n+1)));
end

figure(1); plot(tvals,yvals);

maxerr = max(abs(yvals - tvals.^2));

fprintf('\n nsteps = %d, maxerr = %g \n',nsteps, maxerr);

%% Demonstrates MATLAB's ode45 on the Gear IVP. 

y_0 = 0;
t_0 = 0; 
t_f = 1;

fun = @(t,y)-1000*(y-t^2) + 2*t;

% Set the error tolerances to 1e-5. 
myopts = odeset('RelTol',1e-5,'AbsTol',1e-5);
% Call ode45.
[tvals,yvals] = ode45(fun,[t_0 t_f],y_0,myopts);

maxerr = max(abs(yvals - tvals.^2));
nsteps = size(tvals,1);

fprintf('\n nsteps = %d, no. f-evals = %d, maxerr = %g \n',...
    nsteps, 6*nsteps, maxerr);

figure(1); plot(tvals,yvals);
