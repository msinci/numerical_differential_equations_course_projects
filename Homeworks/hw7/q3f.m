function [ vectorY ] = q3f( t, x )
    vectorY = zeros(2,1);
    vectorY(1) = x(2);
    vectorY(2) = -x(1)-100*x(2);
end

