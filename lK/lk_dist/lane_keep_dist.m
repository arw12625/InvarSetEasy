
%Load system variables
load('A', 'A');
load('B', 'B');
load('E', 'E');
load('U', 'U');
load('X', 'X');
load('W', 'W');
C = E;

%%

%Omega = diag([1,0.1,1,0.1]) * Polyhedron.unitBox(size(A,1));
%Omega = computeOuterApproxInvariant(X,U,10,A,B);
Omega = X;

Nrange = 15:22;

lkTime = zeros(size(max(Nrange),1),1);
betas = zeros(size(max(Nrange),1),1);
for N = Nrange
tic;
[beta, diagnos, sH, sF] = controlInvariantDistGrowthLP(Omega, X, U, W, N, A, B, C);
%[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B);
betas(N) = beta;
lkTime(N) = toc;
end
lkTime = lkTime(Nrange);
betas = betas(Nrange);

betas


%%

s0 = 1 / beta * Omega;
s = s0;
ver = [];

% The invariant set is the convex hull of the backwards step sets
invar = Polyhedron(ver);

for i = 1:N
    s = polyPreDist(A,B,C,X,U,W,s);
    volume(s)
    ver = [ver; s.V];
end

invar = Polyhedron(ver);

verts = s0.V;
satisfied = zeros(size(verts,1), 1);
for i = 1:size(verts, 1)
    satisfied(i) = s.contains(verts(i,:)'); 
end
sum(satisfied) == size(satisfied, 1)

nomVol = volume(invar)

%%

t = X;
iter = 5;

for i = 1:iter
    t = t & polyPre(A,B,X,U,t);
end

overVol = volume(t)