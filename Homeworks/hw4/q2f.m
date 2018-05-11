function [ out ] = q2f(x,y,z)
    a = 10;
    b = 8/3;
    r = 28;

    out = [-a*(y(1)-y(2)); 
        y(1)*(r-y(3))-y(2); 
        -b*y(3)+y(1).*y(2) ] ;
end

