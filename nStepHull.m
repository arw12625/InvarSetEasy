function [invar, S] = nStepHull(Omega, X, U, N, A, B, alpha)
%   Omega - initial set, seed
%   X - State space / constraints
%   U - Input space
%   N - number of steps
%   A - non-singular state-transition matrix
%   B - input matrix
%   alpha - scaler found by growth LP
%
%

s = alpha * Omega;
S = Polyhedron(N);
ver = [];
for i = 1:N
    s = polyPre(A,B,X,U,s);
    S(i) = s;
    ver = [ver; s.V];
end

% The invariant set is the convex hull of the backwards step sets
invar = Polyhedron(ver);

end

