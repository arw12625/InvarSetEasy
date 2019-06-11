function [diagnostics] = testEasyDist(Omega, X, U, W, N, A, B, C)
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

yalmip('clear')

%Linear inequality representation of polytopes X,U,Omega
F = X.A;
f = X.b;
G = U.A;
g = U.b;
H = Omega.A;
h = Omega.b;
E = W.A;
e = W.b;

%Dimensions of system and constraints
n = size(A,2);
m = size(B,2);
l = size(C,2);

%Decision variables for LP
T = sdpvar(size(H,1) + size(F,1) * (N+1) + size(G,1) * N + size(E,1) * N,  size(H,1)+ N* size(E,1), 'full');
Kx = sdpvar(N * m, n, 'full');

Kw = sdpvar(N*m, N*l);
for i = 1:N
    Kw((i-1)*m + (1:(N-i+1)*m),(i-1)*l+(1:l)) = zeros((N-i+1)*m,l);
end

targetInputCon = cell(1,N);
targetDistCon = cell(1,N);
stateStateCon = cell(N+1,1);
stateInputCon = cell(N+1,N);
stateDistCon = cell(N+1,N);

gbar = cell(1+(N+1)+N+N,1);

gbar{1,1} = h;
gbar{1+N+1,1} = f;
stateStateCon{N+1,1} = F;
for i = 1:N
    targetInputCon{1,i} = H * A^(i-1) * B;
    targetDistCon{1,i} = H * A^(i-1) * C;
    stateStateCon{i,1} = F * A^(N-i+1);
    for j = i:N
        stateInputCon{j-i+1,j} = F * A^(i-1) * B;
        stateDistCon{j-i+1,j} = F * A^(i-1) * C;
        if i > 1
            stateInputCon{j,j-i+1} = zeros(size(F,1),m);
            stateDistCon{j,j-i+1} = zeros(size(F,1),l);
        end
        stateInputCon{N+1,i} = zeros(size(F,1),m);
        stateDistCon{N+1,i} = zeros(size(F,1),m);
    end
    gbar{1+i,1} = f;
    gbar{1+(N+1)+i,1} = g;
    gbar{1+(N+1)+N+i,1} = e;
end

Gbarx = cell(4,1);
Gbaru = cell(4,1);
Gbarw = cell(4,1);

Gbarx{1,1} = H * A^N;
Gbaru{1,1} = cell2mat(targetInputCon);
Gbarw{1,1} = cell2mat(targetDistCon);

Gbarx{2,1} = cell2mat(stateStateCon);
Gbaru{2,1} = cell2mat(stateInputCon);
Gbarw{2,1} = cell2mat(stateDistCon);

Gbarx{3,1} = zeros(size(G,1) * N,n);
Gbaru{3,1} = kron(eye(N), G);
Gbarw{3,1} = zeros(size(G,1) * N, l * N);

Gbarx{4,1} = zeros(size(E,1) * N,n);
Gbaru{4,1} = zeros(size(E,1) * N, m * N);
Gbarw{4,1} = kron(eye(N), E);

Gbarx = cell2mat(Gbarx);
Gbaru = cell2mat(Gbaru);
Gbarw = cell2mat(Gbarw);

Hbar = blkdiag(H,kron(eye(N),E));
hbar = [h;repmat(e,N,1)];

gbar = cell2mat(gbar);

Constraints = [
    T*Hbar == [Gbarx + Gbaru * Kx, Gbaru * Kw + Gbarw];
    T*hbar <= gbar;
    T >= 0;
];
Options = sdpsettings('solver', 'gurobi', 'verbose', 1);


diagnostics = optimize(Constraints, [], Options);

tmp = value(T*hbar - gbar);
tmp(size(H,1)+N*size(F,1) + (1:size(F,1)))
value(Kw)

check(Constraints)

end

