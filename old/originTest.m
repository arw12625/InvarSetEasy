U = 1 * Polyhedron.unitBox(2); %Input Set

verts = [2,1;2,4;-2,4;-2,1];
Omega = Polyhedron('V',verts);
verts = [4,0.5;4,8;-4,8;-4,0.5];
X = Polyhedron('V',verts);

theta = pi / 8;
A = 1.1 * [cos(theta), sin(theta); -sin(theta), cos(theta)];
B = eye(2);

N = 8; %Number of backwards steps to use

[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);
beta


%%

S = Polyhedron(N);

[invar, S] = nStepHull(Omega, X, U, N, A, B, 1 / beta);
figure();
hold on;
plot(invar);

for i=1:N
    plot(S(i));
end

plot(1 / beta * Omega);
alpha(0.1);


