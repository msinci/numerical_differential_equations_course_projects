function Fpr_val = fpr_Bratu(u)
%
% Evaluates F'(u) for the Bratu problem 
%
%   (Laplacian)u + lambda*exp*u = 0 in D = [0,1]x[0,1] 
%                             u = 0 on the boundary of D
% 
% The discretization is by finite differences on an mxm regular 
% grid of interior points in [0,1]x[0,1]. 
%
% This routine should not be called until after a call to 
% f_Bratu.m, which sets the Laplacian L_GL and the parameter 
% lambda_GL and also checks to see that n = m^2 for an integer m. 
% 

global L_GL; 
global lambda_GL;

% L_GL is evaluated and passed in from f_Bratu.
% lambda_GL is passed in from newton_bratu_demo.

n = size(u,1);

Fpr_val = L_GL + lambda_GL*spdiags(exp(u),0,n,n);

