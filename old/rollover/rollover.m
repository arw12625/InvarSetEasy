
%Load system variables
load('A', 'A');
load('B', 'B');
load('U', 'U');
load('X', 'X');
n = size(A,1);

Omega = Polyhedron.unitBox(n);

N = 40;

[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B)

%% Volume estimation using random sampling

membershipTest = nStepHullMembership(Omega, X, U, N, A, B, 1 / beta);
9
innerVolEst = uniformSamplingVolEst(membershipTest, X, 10);

innerVolEst


%%

s0 = 1 / beta * Omega;
s = s0;
ver = [];

'Back comp'

for i = 1:N
    i
    s = polyPre(A,B,X,U,s);
    ver = [ver; s.V];
end

'Invar Comp'

invar = Polyhedron(ver);

'Volume Comp'

nomVol = volume(invar)

%%

t = X;
iter = 5;

for i = 1:iter
    t = t & polyPre(A,B,X,U,t);
end

overVol = volume(t)