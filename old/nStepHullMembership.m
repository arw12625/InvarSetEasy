function [membershipTest] = nStepHullMembership(Omega, X, U, N, A, B, alpha)
%   Omega - initial set, seed
%   X - State space / constraints
%   U - Input space
%   N - number of steps
%   A - non-singular state-transition matrix
%   B - input matrix
%   alpha - scaler found by growth LP
%
%

n = size(A,1);
x = sdpvar(n,1,'full');
membershipConstraints = nStepHullConstraints(Omega, X, U, N, A, B, alpha, x);

options = sdpsettings('solver', 'gurobi', 'verbose', 0);

interFunction = @(diagnostics) 1 - diagnostics.problem;
membershipTest = @(xp)  interFunction(optimize([membershipConstraints; x == xp], [], options));

end

