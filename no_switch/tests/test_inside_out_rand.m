%% System definitions
n=8;
m=4;
A = randn(n);
B = randn(n,m);
E = randn(n,n)/n^2;
C = eye(n);
G = zeros(n,0);

sys = LTISys(A,B,E,C,G);

% The horizon considered
T = 8;

% System constraints as polyhedra
X = Polyhedron.unitBox(n);
U = Polyhedron.unitBox(m);
W = 0.5* Polyhedron.unitBox(n);
V = Polyhedron.fullSpace(0);

%% Invariant set computations

% The seed set used for computation
Omega = 0.3 * Polyhedron.unitBox(n);
%Omega = 0.001 * (Xsafe & Polyhedron.unitBox(n));
%Omega = ExamplePoly.poly3d_sin('d',n,'n',5);
%Omega = ExamplePoly.randHrep('d',n,'n',30);
Omega.minHRep;
P = Omega.A;
q0 = zeros(size(P,1),1);
%q0 = q1s;
%q0 = qqq;

uset = Polyhedron.fullSpace(n) * polyhedronPower(W,T) * polyhedronPower(V,T);
aset = polyhedronPower(U,T) * polyhedronPower(X,T) * X;
%uset = Omega * polyhedronPower(W,T) * polyhedronPower(V,T);
%aset = polyhedronPower(U,T) * polyhedronPower(X,T) * Omega;

cf = ControlFamily(sys, T, 1);

options = sdpsettings('verbose', 1, 'solver', 'gurobi'); % options for the LP solver

%c = rand(size(P,1),1)+0.05;
c = ones(size(P,1),1);

[q, converged, infeasible] = insideOutRecurrenceOpt(sys,cf,uset,aset, P, q0, c, 15, 0.01, options);

converged
infeasible

%diagnostics = computeController(sys,cf,uset,aset, options);
%diagnostics.problem