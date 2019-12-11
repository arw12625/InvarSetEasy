function [d] = computeDisturbanceOffsets(plsys, Omega)
%computeDisturbanceOffsets Computes the halfspace offsets that represent
%the effect of the disturbance upon Pre operation
%   plsys - polytopic linear system
%   Omega - polytope to compute pre of
%
%   d - disturbance offset

N = 1;

F = Omega.A;
f = Omega.b;
E = plsys.W.A;
e = plsys.W.b;

if N <= 0
    d = zeros(size(F,1), 1);
    return;
end

if sum(sum(abs(-F * plsys.A^(N-1) * plsys.E))) == 0
    d = zeros(size(F,1), 1);
    return;
end

w = sdpvar(size(plsys.E,2),1,'full');

objective = -F * plsys.A^(N-1) * plsys.E * w;

constraints = [ E * w <= e ];

options = sdpsettings('solver', 'gurobi', 'verbose', 0);

diagnostics = optimize(constraints, objective, options);

d = zeros(size(F,1),1);

for i = 1:size(F,1)
    selectsolution(i);
    d(i) = value(objective(i));
end

end

