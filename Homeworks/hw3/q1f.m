function [ out ] = q1f(t, y)

D = 0.1;
out = [y(2); -y(1)-D*y(2)];

end

