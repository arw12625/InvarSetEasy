%% System definitions

% This script demonstrates various invariant set computations for a 
% discrete time double integrator using this library. 

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + w(t)
% y(t) = C x(t) + v(t)
A = [1,1;0,1];
B = [0;1];
E = [1;0.1];
C = [1,0];
G = zeros(1,0);

sys = LTISys(A,B,E,C,G);

% The horizon considered
T = 6;

% System constraints as polyhedra
Xsafe = (1 * Polyhedron.unitBox(1)) * (1 * Polyhedron.unitBox(1));
Usafe = Polyhedron.unitBox(1);
W = 0.125 * Polyhedron.unitBox(1);
V = Polyhedron.fullSpace(0);

%% Invariant set computations

% The seed set used for computation
Omega = 0.2 * Polyhedron.unitBox(2);

uset = Xsafe * polyhedronPower(W,T) * polyhedronPower(V,T);
aset = polyhedronPower(Usafe,T) * polyhedronPower(Xsafe,T) * Xsafe;

cf = ControlFamily(sys, T, 0);

options = sdpsettings('verbose', 0, 'solver', 'gurobi'); % options for the LP solver

c = ones(size(Omega.b));
[q, converged, infeasible] = insideOutRecurrenceOpt(sys,cf,uset,aset, Omega.A, Omega.b, c, 5, 0.001, options);
converged
infeasible