function [d] = computeDisturbanceOffsets(X, W, N, A, C)
%computeDisturbanceOffsets Computes the halfspace offsets that represent
%the effect of the disturbance upon Pre operation
%   X - State space / constraints
%   W - Disturbance space
%   N - number of steps
%   A - non-singular state-transition matrix
%   C - disturbance matrix
%
%   d - disturbance offset

F = X.A;
f = X.b;
E = W.A;
e = W.b;

if N <= 0
    d = zeros(size(F,1), 1);
    return;
end

w = sdpvar(size(C,2),1,'full');

objective = -F * A^(N-1) * C * w;

constraints = [ E * w <= e ];

options = sdpsettings('solver', 'gurobi', 'verbose', 0);

diagnostics = optimize(constraints, objective, options);

d = zeros(size(F,1),1);

for i = 1:size(F,1)
    selectsolution(i);
    d(i) = value(objective(i));
end

end

