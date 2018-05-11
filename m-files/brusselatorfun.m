function fval = brusselatorfun(t,y)
%
% Evaluates the right-hand side of the Brusselator ODE 
%
%  (y_1)' = A + y_1^2*y_2 - (B+1)*y_1
%
%  (y_2)' = B*y_1 - y_1^2*y_2
%
% From E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary 
% Differential Equations I", 2nd ed., Springer Series in Computational 
% Mathematics, Vol. 8, 1993. 
%

A = 1; B = 3;

fval = [A + (y(1,1)^2)*y(2,1) - (B+1)*y(1,1); ... 
        B*y(1,1) - (y(1,1)^2)*y(2,1)]; 
