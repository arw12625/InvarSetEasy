function [diagnostics] = testStateControlInvariance(Omega, X, U, N, A, B)
%testStateControlInvariance Determine if the N backwards step set of 
%(Omega,X,U) is control invariant with state constraints
%   Omega - initial set, seed
%   X - State space / constraints
%   U - Input space
%   N - number of steps
%   A - non-singular state-transition matrix
%   B - input matrix
%
%   This is taken from Theorem 2 in Computing control invariant sets in high dimension is easy
%   Omega,X,U are convex polytopes each containing the origin
%       They should have halfspace representations in MPT
%   This function requires YALMIP
%
%   For now ignore equality constraints in polytopes but this can be added
%   Also this function might apply to singular A matrices without modification
%

%Linear inequality representation of polytopes X,U,Omega
F = X.A;
f = X.b;
G = U.A;
g = U.b;
H = Omega.A;
h = Omega.b;

%Dimensions of system and constraints
n = size(A,1);
m = size(B,2);
nf = size(F,1);
ng = size(G,1);
nh = size(H,1);

%Dimensions of LP
ngt = nh + N * ng + (N+1) * nf ;
nt = n + N * m;

%Decision variables for LP
T = sdpvar(ngt,  nh, 'full');
M = sdpvar(nt, nt, 'full');

%Constraints for LP
Gt = zeros(ngt, nt);
gt = zeros(ngt,1);
Ht = [H, zeros(nh, nt-n)];

Gt(1:nh, 1:n) = H * A^N;
gt(1:nh) = h;
for i = 1:N
    Gt(1:nh, n + (i-1)*m + (1:m)) = H * A^(i-1) * B;
    Gt(nh + (i-1)*ng + (1:ng), n + (i-1)*m + (1:m)) = G;
    Gt(nh + nh * N + (i-1)*nf + (1:nf), (1:n)) = F * A^(N-i+1);
    for j = i:N
        Gt(nh + ng * N + (i-1)*nf + (1:nf), n + (j - 1) * m + (1:m)) = F * A^(j-i) * B;
    end
    gt(nh + (i-1)*ng + (1:ng)) = g;
    gt(nh + N * ng + (i-1)*nf + (1:nf)) = f;
end
Gt((ngt - nf + 1):ngt, 1:n) = F;
gt((ngt - nf + 1):ngt) = f;

identMat = [eye(n), zeros(n, nt - n)];

Constraints = [
    T*Ht == Gt * M;
    T*h <= gt;
    T >= 0;
    identMat * M == identMat;
];

%   Run feasibility LP
diagnostics = optimize(Constraints);


end

