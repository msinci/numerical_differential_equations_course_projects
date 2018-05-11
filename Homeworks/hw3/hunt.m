function [ out ] = hunt(t, y)

alpha = 0.2;
out = [2*y(1)-alpha.*y(1).*y(2); -y(2)+alpha.*y(1).*y(2)];

end

