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
C = [1,0];

sys = LTISys(A,B,C);

% The horizon considered
T = 6;

% System constraints as polyhedra
Xsafe = (1 * Polyhedron.unitBox(1)) * (1 * Polyhedron.unitBox(1));
Usafe = Polyhedron.unitBox(1);
Xterm = Xsafe;
W = [0; 0.125] * Polyhedron.unitBox(1);
V = 0.125 * Polyhedron.unitBox(1);

%% Invariant set computations

% The seed set used for computation
Omega = (0.2* Polyhedron.unitBox(1)) * (0.2 * Polyhedron.unitBox(1));

uset = Omega * polyhedronPower(W,T) * polyhedronPower(V,T);
aset = polyhedronPower(Xsafe,T) * Xterm * polyhedronPower(Usafe,T);

cf = ControlFamily(sys, T);

options = sdpsettings('verbose', 1, 'solver', ''); % options for the LP solver

diag = computeController(sys,cf,uset,aset,options);