function B = biharmonic_op(m)
%
% Returns the matrix for the biharmonic operator on a regular mxm grid 
% of interior points, assuming the standard treatment of the normal 
% derivative boundary conditions, i.e., reflecting the first interior 
% gridline across the boundary on each side. 
%
% The input is  
%  m        = size of each grid line. 
%
I = speye(m,m);
E1 = sparse(2:m,1:m-1,1,m,m);
E2 = sparse(3:m,1:m-2,1,m,m);
T = 20*I - 8*(E1 + E1') + E2 + E2'; T(1,1) = 21; T(m,m) = 21; 
S = -8*I + 2*(E1 + E1'); 
B = kron(I,T) + kron(E1 + E1',S) + kron(E2 + E2',I);
B(1:m,1:m) = B(1:m,1:m) + I; 
B(m^2-m+1:m^2,m^2-m+1:m^2) = B(1:m,1:m); 
