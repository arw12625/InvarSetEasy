% This ran for a random 20 dimensional system with 20 steps in 30 minutes


n = 20;
A = (rand(n) - 0.5) * 1.5;
B = eye(n);
%C = eye(n);

X = 100 * Polyhedron.unitBox(n);
U = 10 * Polyhedron.unitBox(n);
%W = 0.1 * Polyhedron.unitBox(n);
Omega = Polyhedron.unitBox(n);

N = 20;

tic;
[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);
scaleTime = toc;

scaleTime

%[beta2, diag2, ~, ~] = controlInvariantDistGrowthLP(Omega, X, U, W, N, A, B, C)