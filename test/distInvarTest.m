%% Compute invariant set scaling

n = 2; %System Dimension
Omega = 2 * Polyhedron.unitBox(n);
X = 5 * Polyhedron.unitBox(n);
U = 0.6 * Polyhedron.unitBox(n);
W = 0.05 * Polyhedron.unitBox(n);

theta = pi / 4;
A = 1.1 * [cos(theta), sin(theta); -sin(theta), cos(theta)];
B = eye(n);
E = eye(n);

N = 2;

plsys = PolyLinSys(A,X,B,U,E,W);

yalmipOptions = sdpsettings('verbose', 1);

[scale, isInvar, diagnostics, Kx,Kw,u0] = growEasyControlInvariance(Omega, plsys, N, 1, 0, yalmipOptions);
beta = 1 / scale;
%[beta, diag, sF, sH] = controlInvariantDistGrowthLPMOD(Omega, X, U, W, N, A, B, C);

alphaS = 1 / beta;

%% Compute and plot the scaled backwards step sets explicitly

figure()
hold on;
s = alphaS * Omega;
S = Polyhedron(N);
ver = [];
for i = 1:N
    s = polyLinPre(plsys,s);
    S(i) = s;
    plot(s);
    ver = [ver; s.V];
end

% The invariant set is the convex hull of the backwards step sets
invar = Polyhedron(ver);
plot(invar, 'Color', 'b');

%Lower opacity to visualize overlap
alpha(0.1);

plot(alphaS * Omega, 'Color', 'g');
