function seesol(u,m)
%
% Forms a surface plot for the solution of a Dirichlet problem on the unit 
% square in R^2 discretized on an evenly spaced mesh of mxm interior points. 
%
% The inputs are 
%   u       = solution values at mxm interior mesh points expressed as a 
%             vector of length m^2. 
%   m       = number of interior mesh points per side. 
% 
h = 1/(m+1);
x = 0:h:1;
y = x;
ubar = reshape(u,m,m);
z = zeros(m,1);
ubar = [z,ubar,z];
z = zeros(m+2,1);
ubar = [z,ubar',z];
surf(x,y,ubar);
%mesh(x,y,ubar);
view(5,45)
