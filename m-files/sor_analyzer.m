function sor_analyzer(A)
%
% First evaluates and prints the spectral radius of the Gauss-Seidel 
% iteration matrix. Then evaluates and plots the spectral radius of 
% the SOR iteration matrix as a function of omega over the interval 
% (0,2) and prints the optimal omega and spectral radius.
%
% The input is 
%  A       = matrix
%
% M-file required:
%   spec_rad.m               (computes the spectral radius)
%

L=tril(A,-1);
U=triu(A,1);
D=diag(diag(A)); 

% Evaluate and print the spectral radius of the Gauss-Seidel 
% iteration matrix. 
TGS = (D+L)\U; rho_GS = spec_rad(TGS); 
fprintf('\n For Gauss-Seidel, spectral radius = %f \n',rho_GS); 

% Evaluate and plot the spectral radius of the SOR iteration 
% matrix as a function of omega
h = 101\2;
w = h:h:2-h; rho = zeros(size(w)); 
for i = 1:length(w)
  TSOR = (D+w(i)*L)\( (1-w(i))*D - w(i)*U );
  rho(i) = spec_rad(TSOR);
end
[rho1,ind] = sort(rho); 
rho_star = rho1(1);w_star = w(ind(1));
fprintf('\n For SOR, optimal relaxation parameter = %f,  optimal spectral radius = %f \n',w_star,rho_star);
plot(w,rho);grid on;
%axis([0 2 0 1]);
drawnow;
