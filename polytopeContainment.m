function [diagnostics] = polytopeContainment(Gamma, Delta)
%polytopeContainment tests if polytope Gamma is contained in polytope Delta
%   Gamma,Delta - Polytopes for testing
%       Polytopes should be specified in halfspace form

%   For now ignore equality constraints in polytopes
H = Gamma.A;
h = Gamma.b;
G = Delta.A;
g = Delta.b;

n = size(G,2);
q = size(G,1);
p = size(H,1);

T = sdpvar(q,p, 'full');

Constraints = [
    T*H == G;
    T*h <= g;
    T >= 0;
];

diagnostics = optimize(Constraints);


end

