function [xf] = stepMap(steps, A, B, x0, u)
%stepMap computes the trajectory of a point in a linear system with input
%   steps - number of steps
%   A - state matrix
%   B - input matrix
%   x0 - initial point
%   u - input

xf = A^steps * x0;
for step = 0:(steps - 1)
    xf = xf + A^(steps - step - 1) * B * u(:, steps - step);
end

end

