function w = adi_pre(v)
%
% Implements a specified number of sweeps of 
% alternating direction-implicit (ADI) iteration. 
%
% Input: 
%   v      = vector
%
% Output:
%   x      = result of applying the ADI sweeps to v. 
%
% The "horizontal" and "vertical" operators H_GL and V_GL are set
% in the calling program and passed in as global variables. 

global H_GL V_GL 

% NOTE: The following are naive choices of rho. There may be room
% for improvement. 
rho = .05; 
%rho = 1/sqrt(size(H_GL,1));

rhoI = rho*speye(size(H_GL)); 
w = zeros(size(v)); 
nsweeps = 1; 
for i = 1:nsweeps
   w = (H_GL + rhoI)\((rhoI - V_GL)*w + v); 
   w = (V_GL + rhoI)\((rhoI - H_GL)*w + v);
end
