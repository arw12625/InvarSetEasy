function [beta, diagnostics, sF, sH] = controlInvariantDistGrowthLP(Omega, X, U, W, N, A, B, C)
%controlInvariantDistGrowthLP finds a scaling of Omega so that the
%N backward step set computed is control invariant with disturbance.
%   Omega - initial set, seed
%   X - State space / constraints
%   U - Input space
%   W - Disturbance space
%   N - number of steps
%   A - non-singular state-transition matrix
%   B - input matrix
%   C - disturbance matrix
%
%   This is motivated by Theorem 2 in 'Computing control invariant sets in high dimension is easy'
%   I need to finish a formal proof that this method with disturbance is correct
%   Omega,X,U,W are convex polytopes each containing the origin
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
gamma = sdpvar(1);

% gamma is the reciprocal of the scaling applied to Omega
% so to obtain the largest invariant set, gamma must be minimized
Objective = gamma;

% Compute effect of disturbance through offsets
dH = zeros(size(H,1), N+1);
dF = zeros(size(F,1), N+1);
sF = zeros(size(F,1), N+1);
for i = 1:(N+1)
    dH(:,i) = computeDisturbanceOffsets(Omega,W,i-1,A,C);
    dF(:,i) = computeDisturbanceOffsets(X,W,i-1,A,C);
    if i > 1
        sF(:,i) = sF(:,i-1) + dF(:,i);
    elseif i == 1
        sF(:,i) = dF(:,i);
    end
end
sH = sum(dH,2);

%Constraints for LP
Gt = zeros(ngt, nt);
gh = zeros(ngt,1);
gtild = zeros(ngt,1);
gtild(1:nh) = h;
gh(1:nh) =  sH;
Ht = [H, zeros(nh, nt-n)];

Gt(1:nh, 1:n) = H * A^N;
for i = 1:N
    Gt(1:nh, n + (i-1)*m + (1:m)) = H * A^(i-1) * B;
    Gt(nh + (i-1)*ng + (1:ng), n + (i-1)*m + (1:m)) = G;
    Gt(nh + ng * N + (i-1)*nf + (1:nf), (1:n)) = F * A^(N-i+1);
    for j = i:N
        Gt(nh + ng * N + (i-1)*nf + (1:nf), n + (j - 1) * m + (1:m)) = F * A^(j-i) * B;
    end
    gh(nh + (i-1)*ng + (1:ng)) = g;
    gh(nh + N * ng + (i-1)*nf + (1:nf)) = f + sF(:,N - i + 2);
end
Gt((ngt - nf + 1):ngt, 1:n) = F;
gh((ngt - nf + 1):ngt) = f + sF(:,1);

identMat = [eye(n), zeros(n, nt - n)];

Constraints = [
    T*Ht == Gt * M;
    T*(h ) <= gamma * gh + gtild;
    T >= 0;
    gamma >= 0;
    identMat * M == identMat;
];

Options = sdpsettings('solver', 'gurobi', 'verbose', 1);

diagnostics = optimize(Constraints, Objective, Options);

beta = value(gamma);
end

