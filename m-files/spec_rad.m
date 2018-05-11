function rho = spec_rad(T)
%
% Computes the spectral radius of the matrix T. 
%
% The input is 
%  T       = square matrix (can be full or sparse) 
%
% The output is 
%  rho     = spectral radius of T
%
if issparse(T)
  options.disp = 0;
  rho = abs(eigs(T,1,'lm',options));
else
  rho = max(abs(eig(T)));
end
