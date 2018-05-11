function mol_demo_1D(m,tInt,funIC)
%
% Demonstrates method-of-lines (MoL) solution of the PDE IBVP 
%
%     u_t = u_xx,           0 < x < 1, 0 < t
%     u(0,t) = u(1,t) = 0,  0 < t, 
%     u(x,0) = f(x),        0 < x < 1.
%
% Alternatively (by uncommenting appropriate lines) demonstrates MoL
% solution of the nonlinear IBVP 
%
%     u_t = u_xx + u^2,           0 < x < 1, 0 < t
%     u(0,t) = u(1,t) = 0,  0 < t, 
%     u(x,0) = f(x),        0 < x < 1.
%
% Dicretization with respect to x is done on an evenly 
% spaced mesh of m interior points in [0,1]. Both discretized problems
% become very stiff very quickly as m grows. 
%
% Inputs:
%   m       = number of interior mesh points in [0,1].
%   tInt    = row vector [t_0 t_f] of initial and final times.
%   funIC   = function handle (@funIC) for the initial-value 
%	          function f. This should accept an argument x 
%	          and return f(x). 
%

% Set up the step and the mesh points. 
h = (m+1)\1;
x = (h:h:1-h); 

% Set up the initial values.
u_0 = zeros(m,1);
for i = 1:m
  u_0(i,1) = funIC(x(1,i));
end

figure(gcf);plot([0,x,1],[0,u_0',0]);
%axis([0 1 min(min(u_0)) 1.2*max(max(u_0))]);
axis([0 1 0 1.2*max(max(u_0))]);
title(['Initial Function, m = ' num2str(m)]);
disp('Strike any key to continue.');
pause;

% Form the sparse discrete 2nd-derivative matrix. 
I = speye(m,m);
E = sparse(2:m,1:m-1,1,m,m);
L =  (h^2)\(-2*I + E + E');

% Define the right-hand side of the ODE and the appropriate options
% structure. 
% Uncomment the next two lines for u_t = u_xx:
% rhsfun = @(t,u)L*u;
% options = odeset('Stats','on','Jacobian',L);
% Uncomment the next two lines for u_t = u_xx + u^2:
rhsfun = @(t,u)L*u + u.^2;
options = odeset('Stats','on','Jacobian',@(t,u)L+2*spdiags(u,0,m,m));

% Call the ODE solver. 
% Uncomment the next line for ode45 (very inefficient).  
% [T,Y] = ode45(rhsfun,tInt,u_0,options); method = 'ode45';
% Uncomment the next line for ode15s (much more efficient). 
 [T,Y] = ode15s(rhsfun,tInt,u_0,options); method = 'ode15s';

figure(gcf);

for i = 1:size(T,1)
  plot([0,x,1],[0,Y(i,:),0]);axis([0 1 min(min(Y)) 1.2*max(max(Y))]);
  title([method, ' time = ', num2str(T(i,1))]);drawnow;pause(.1);
end
