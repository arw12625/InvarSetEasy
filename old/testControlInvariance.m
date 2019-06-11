function [diagnostics] = testControlInvariance(Omega, U, N, A, B)
%testControlInvariance Test if the hull of the N backwards step set of
%(Omega,X,U) is control invariant
%   Omega - initial set, seed
%   U - Input space
%   N - number of steps
%   A - non-singular state-transition matrix
%   B - input matrix
%
%   This is taken from Theorem 1 in 'Computing control invariant sets is easy'
%   Omega,X,U are convex polytopes each containing the origin
%       They should have halfspace representations in MPT

%   For now ignore equality constraints in polytopes
G = U.A;
g = U.b;
H = Omega.A;
h = Omega.b;

n = size(A,1);
m = size(B,2);
ng = size(G,1);
nh = size(H,1);

ngt = nh + N * ng;
nht = nh + 2 * N * m;
nt = n + N * m;

Gt = zeros(ngt, nt);
gt = zeros(ngt,1);
Ht = zeros(nht, nt);
ht = zeros(nht, 1);

Gt(1:nh, 1:n) = H * A^N;
gt(1:nh) = h;
Ht(1:nh, 1:n) = H;
ht(1:nh) = h;

for i = 1:N
    Gt(1:nh, n + (i-1)*m + (1:m)) = H * A^(i-1) * B;
    Gt(nh + (i-1)*ng + (1:ng), n + (i-1)*m + (1:m)) = G;
    gt(nh + (i-1)*ng + (1:ng)) = g;
    Ht(nh + (i-1)*m + (1:2*m), n + (i-1)*m + (1:m)) = [eye(m); -eye(m)];
end

T = sdpvar(ngt,  nht, 'full');
M = sdpvar(nt, nt, 'full');
identMat = [eye(n), zeros(n, nt - n)];

Constraints = [
    T*Ht == Gt * M;
    T*ht <= gt;
    T >= 0;
    identMat * M == identMat;
];

diagnostics = optimize(Constraints);

end

