function [ vectorY ] = q4f( t, x )
    vectorY = zeros(2,1);
    vectorY(1) = 77.27*( x(2) + x(1)*(1- 8.375*10^-6 *x(1) -x(2)));
    vectorY(2) = (x(3)-(1+x(1)*x(2)))/77.27;
    vectorY(3) = 0.161*(x(1)-x(3));
end

