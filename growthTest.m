%% Compute invariant set scaling

n = 2; %System Dimension
Omega = Polyhedron.unitBox(n); %
X = 5 * Polyhedron.unitBox(n);
U = 0.6 * Polyhedron.unitBox(n);

theta = pi / 4;
A = 1.2 * [cos(theta), sin(theta); -sin(theta), cos(theta)];
B = eye(n);

N = 6;

figure()
hold on;

[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);
beta
diag

alphaS = 1 / beta;

%% Compute and plot the scaled backwards step sets explicitly

s = alphaS * Omega;
S = Polyhedron(N);
ver = [];
for i = 1:N
    s = polyPre(A,B,X,U,s);
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


%% Test random points for membership in invariant set using given constraints

x = sdpvar(n, 1, 'full');
Constraints =  nStepHullConstraints(Omega, X, U, N, A, B, alphaS, x);

testNum = 100;
cols = zeros(testNum,1);
testPoints = 4 * rand(2,testNum) - 2;

for i = 1:testNum

testConstraints = [Constraints; x == testPoints(:,i)];

cols(i) = 3;
diag = optimize(testConstraints);
if diag.problem == 0
    cols(i) = 1;
end

end
hold on;
scatter(testPoints(1,:), testPoints(2,:), [], cols);

% Note - There appears to be an error or inaccuracy in this method.
% Some points contained in the invariant set near the boundary are infeasible in the LP