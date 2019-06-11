
%Load system variables
load('A', 'A');
load('B', 'B');
load('U', 'U');
load('X', 'X');

Omega = X;

N = 15;

tic;
[beta, diag] = controlInvariantGrowthLP(Omega, X, U, N, A, B)
lkTime = toc;


%For 15 steps, growth took 0.59 seconds and the corresponding invariant set
%had at least 99.72% of the volume of the maximal set


%%

s0 = 1 / beta * Omega;
s = s0;
ver = [];

% The invariant set is the convex hull of the backwards step sets
invar = Polyhedron(ver);

for i = 1:N
    s = polyPre(A,B,X,U,s);
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