function fval = oregonatorfun(t,y)
%
% Evaluates the right-hand side of the oregonator ODE 
%
%  (y_1)' = 77.27*[y_2 + y_1*(1 - 8.375e-6*y_1 - y_2)]
%
%  (y_2)' = (77.27\1)*[y_3 - (1 + y_1)*y_2]
%
%  (y_3)' = 0.161*(y_1 - y_3) 
%
% From E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary 
% Differential Equations I", 2nd ed., Springer Series in Computational 
% Mathematics, Vol. 8, 1993. 
%

fval = [77.27*(y(2,1) + y(1,1)*(1 - (8.375e-6)*y(1,1) - y(2,1))); ...
        77.27\(y(3,1) - (1 + y(1,1))*y(2,1)); ...
	0.161*(y(1,1) - y(3,1))];
