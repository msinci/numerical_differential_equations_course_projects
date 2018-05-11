%% Runs newtonbt.m on the 2D Bratu problem
%
%  Laplacian*u + lambda*exp(u) = 0 in D = [0,1]x[0,1]
%                            u = 0 on the boundary of D
%
% 
% M-files required: 
%   laplacian.m   (sets up the discretized laplacian)
%   loadrhs.m     (loads the right-hand side)
%   newtonbt.m    (runs Newton's method with backtracking)
%   rhsfun.m      (user-supplied, called by loadrhs.m) 
%   seesol.m      (plots the solution surface)
%

% Set the mesh size
m = 64;

% Set the problem parameter, declaring it global to pass it to the 
% function- and Jacobian-evaluation routines f_Bratu and fpr_Bratu. 
global lambda_GL;
lambda_GL = 2;


% Set the stopping criteria for newtonbt. 
tol_F = 1.e-9;
tol_x = 1.e-9;
itmax = 100;

% First, start from the zero initial guess. 

u0 = zeros(m^2,1);

u = newtonbt(@f_Bratu,@fpr_Bratu,u0,tol_F,tol_x,itmax);

% Plot the solution in figure(1).
figure(1); seesol(u,m);
title(['Bratu solution, lambda = ' num2str(lambda_GL) ...
    ', zero initial guess']);

% Now start from a non-zero initial guess to find the other solution. 

% "Tent" function. 
u0fun = @(x,y)6*(x*(1-x) + y*(1-y));
% "Peaked" function. 
%u0fun = @(x,y)45*x*(1-x)*y*(1-y);

% Load the initial guess into a vector using the loadrhs routine. 
u0 = loadrhs(u0fun,m);

u = newtonbt(@f_Bratu,@fpr_Bratu,u0,tol_F,tol_x,itmax);

% Plot the solution in figure(2).
figure(2); seesol(u,m);
title(['Bratu solution, lambda = ' num2str(lambda_GL)]);

% For comparison, plot the initial guess in figure(3). 
figure(3); seesol(u0,m);
title(['Initial guess in Figure 2, lambda = ' num2str(lambda_GL)]);
