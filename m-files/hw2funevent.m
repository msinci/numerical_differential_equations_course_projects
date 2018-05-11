function [value,isterminal,direction] = hw2funevent(t,y)
% Defines the event for terminating integration of the HW2 IVP 
% when the derivative goes from positive to negative. 
value = y + exp(t)*cos(t); % Detect when the derivative crosses zero. 
isterminal = 1; % Stop when this occurs, but ...  
direction = -1; % only when crossing from positive to negative. 
