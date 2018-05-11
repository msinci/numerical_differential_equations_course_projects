function Jval = brusselatorfunjac(t,y)
%
% Evaluates the Jacobian with respect to y of the Brusselator ODE 
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

Jval = [2*y(1,1)*y(2,1) - (B+1), y(1,1)^2; ... 
        B - 2*y(1,1)*y(2,1), -y(1,1)^2]; 
