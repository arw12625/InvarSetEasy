% Finding a controlled invariant set via a one-shot method

% System matrices, state = [x; dx]
A = [0 1;
     0 0];
B = [0;1];
F = 0 * [0.1;0];
% Bounds
xmax = 1;
umax = 1;

% Time discretization
dt = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adt = expm(A*dt);
Bdt = dt*B;
Fdt = dt*F;

XU = Polyhedron('H', [0 0 1 umax;
				      0 0 -1 umax]);

d = Dyn(Adt, Fdt, Bdt, XU);

S = Polyhedron('A', [eye(2); -eye(2)], 'b', [xmax; xmax; xmax; xmax]);

N = 3;

% Compute attractor defining invariant set
X0 = d.win_always_oneshot(S, N, 0.15);

%%

U = Polyhedron('H', [1, umax; -1, umax]);

[beta, diag] = controlInvariantGrowthLP(X0, S, U, N, Adt, Bdt);

%%

% testStateControlInvariance(1 / (beta * 0.99) * X0, S, U, N, A, B)

%%
X1 = X0;
X2 = 1 / beta * Omega;

for i=1:N
	X1 = intersect(S, d.pre(X1));
	X2 = intersect(S, d.pre(X2));
end

clf;
figure();
hold on;
plot(X1);
plot(X2);
alpha(0.1);