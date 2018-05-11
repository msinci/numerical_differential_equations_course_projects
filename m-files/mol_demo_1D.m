function mol_demo_1D(m,tint,funIC)
%
% Demonstrates method-of-lines solution of the PDE IBVP 
%
%     u_t = u_xx,           0 < x < 1, 0 < t
%
%     u(0,t) = u(1,t) = 0,  0 < t, 
%
%     u(x,0) = f(x),        0 < x < 1.
%
% Dicretization with respect to x is done on an evenly 
% spaced mesh of m interior points in [0,1]. 
%
% Inputs:
%   m       = number of interior mesh points in [0,1].
%   tint    = row vector [t_0 t_f] of initial and final times.
%   funIC   = function handle (@funIC) for the initial-value 
%	      function f. This should accept an argument x 
%	      and return f(x). 
%

% Set up the step and the mesh points. 
h = (m+1)\1;
x = (h:h:1-h); 

% Set up the initial values.
u_0 = zeros(m,1);
for i = 1:m
  u_0(i,1) = feval(funIC,x(1,i));
end

% Form the sparse discrete 2nd-derivative matrix. 
I = speye(m,m);
E = sparse(2:m,1:m-1,1,m,m);
L =  (h^2)\(-2*I + E + E');

% Define the right-hand side of the ODE. 
rhsfun = @(t,u)L*u;

% Set up the ODE solver options.
options = odeset('Stats','on','Jacobian',L);

% Call the ODE solver. 
[T,Y] = ode15s(rhsfun,tint,u_0,options);
%[T,Y] = ode45(rhsfun,tint,u_0,options);

for i = 1:size(T,1)
  plot([0,x,1],[0,Y(i,:),0]);axis([0 1 min(min(Y)) 1.2*max(max(Y))]);
  title(['Time = ',num2str(T(i,1))]);drawnow;pause(.1);
end
