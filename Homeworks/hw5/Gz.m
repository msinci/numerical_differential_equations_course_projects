function [ out ] = Gz( z )
    % implementation of the g(z) to draw the absolute 
    % stability region of the RK4
    out = abs( 1 + z + (1/2)*(z.^2) + (1/6)*(z.^3) + (1/24)*(z.^4)) -1; 
end

