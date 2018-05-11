function [ value, isterminal, direction ] = hw3funevent(t, y)
% define termination conditions
% in our case, if population is below 1
value = 1;
if y(1)<1
    value = 0;
elseif y(2)<1
    value = 0;
end

isterminal = 1; % stop condition
direction = -1; % only when crossing from positive to nefative
