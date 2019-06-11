%% Setup

n = 2; %System Dimension
Omega = Polyhedron.unitBox(n);
X = 5 * Polyhedron.unitBox(n);
U = 0.6 * Polyhedron.unitBox(n);
W = 0.20 * Polyhedron.unitBox(n);

theta = pi / 4;
A = 1.1 * [cos(theta), sin(theta); -sin(theta), cos(theta)];
B = eye(n);
C = eye(n);

N = 5;

plsys = PolyLinSys(A,X,B,U);

%% Compute invariant set with no dist using old method

[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);

%% Compute invariant set scaling


%testEasyDist(1 / beta * Omega,X,U,W,N,A,B,C);
[beta, diagnostics, sF, sH] = controlInvariantDistGrowthLPMOD(Omega,X,U,W,N,A,B,C);


%%

S(N+1) = Polyhedron();
S(1) = 1 / beta * Omega; 
for i = 1:N
    S(i+1) = polyPreDist(A,B,C,X,U,W,S(i));
end

clf
hold on;
plot(S(2:end));
alpha(0.1);
plot(S(1));