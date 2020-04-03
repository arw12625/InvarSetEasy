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
E = [1;0];
C = eye(2);
G = zeros(2,0);

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
Omega = 0.4 * Polyhedron.unitBox(2);

uset = Omega * polyhedronPower(W,T) * polyhedronPower(V,T);
aset = polyhedronPower(Xsafe,T) * Omega * polyhedronPower(Usafe,T);

cf = ControlFamily(sys, T, 1);

options = sdpsettings('verbose', 1, 'solver', ''); % options for the LP solver

diag = computeController(sys,cf,uset,aset,options);

%%

[Ky,uc] = cf.computeOutputFeedbacks();
sysmap = sys.computeSystemMap(T, Ky, uc);
x0 = [-0.0;0.0];
figure(1)
hold on
plot(Xsafe)
plot(Omega)
for i = 1:200
    wcol = 0.25 * rand(T,1) - 0.125;
    vcol = zeros(0,1);
    [x,u,w,v] = sys.simulateClosedLoopRealization(sysmap,x0,wcol,vcol);
    plot(x(1,:),x(2,:))
end
