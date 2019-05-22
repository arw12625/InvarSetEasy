%% Load system variables

load('A', 'A');
load('B', 'B');
load('U', 'U');
load('X', 'X');
n = size(A,1);

Omega = Polyhedron.unitBox(n);

%% Compute scaling for inner approx

N = 100;

tic;
[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);
scaleTime = toc;


%% Compute inner constraints

tic;
innerX = sdpvar(n, 1);
innerConstraints = nStepHullConstraints(Omega, X, U, N, A, B, 1 / beta, innerX);
innerConstraintTime = toc;

%% Compute outer constraints

tic
outerX = sdpvar(n, 1);

numStep = 13;
outerInvarPoly = computeOuterApproxInvariant(X,U,numStep, A, B);
outerConstraints = outerInvarPoly.A * outerX <= outerInvarPoly.b;
outerConstraintTime = toc;

%% Evaluate approximationz

numTests = 20;
origin = zeros(n,1);
[ratios, valid] = evalLinearContainment(innerConstraints,innerX,outerConstraints,outerX,numTests,origin);

ratios
