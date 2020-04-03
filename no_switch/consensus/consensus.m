%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
T = 10; % number of steps
d = 4; % number of UAVs
n = 2 * d; % dimension of state space
m = 1; % dimension of input space
l = d; % dimension of disturbance space
ns = 1; % number of switching modes

Kh = 1;
Kv = 1;

topo = createGridGraph(2,2);
%{
topo(4,5) = 0;
topo(5,4) = 0;
topo(2,3) = 0;
topo(3,2) = 0;
%}

%{
g = digraph(topo{1});
bins = conncomp(g, 'Type', 'weak');
isConnected = all(bins == 1)
%}


Mh = 10;
offset = 0;
Mv = 5;
uH = 2;
wH = 0.1;

%{
% used for 2x7
Mh = 100;
offset = 0;
Mv = 50;
uH = 20;
wH = 0.1;
%}

%{
Mh = 5;
offset = 0;
Mv = 5;
uH = 1;
wH = 0.1;
%}
degrees = sum(topo, 1);
Aol = kron(eye(d), [1,1;0,1]);
Acl = zeros(size(Aol));
for i = 1:d
    if degrees(i) == 0
        continue
    end
    scale = 1 / (degrees(i));
    for j = 1:d
        if j == i
            Acl(2*i, 2*i-1) = -Kh;
            Acl(2*i, 2*i)   = -Kv*(1+scale);
        else
            Acl(2*i, 2*j-1) =  scale * topo(j,i) * Kh;
            Acl(2*i, 2*j)   =  scale * topo(j,i) * Kv;
        end
    end
end

A = Aol + Acl;

%B{i} = [0;1;zeros(n-2, 1)];
%B = [zeros(floor(n/2)-1, 1); 0;1;zeros(ceil(n/2)-1, 1)];
B = repmat([0;1], d,1);
%B{i} = kron(eye(d), [0;1]);

%C = [zeros(2,floor(n/2)-1), eye(2), zeros(2,ceil(n/2)-1)];
%G = zeros(2,0);

C = eye(n);
G = zeros(n,0);

E = kron(eye(d),[0;1]);
%E = zeros(n,0);
%E{i} = zeros(n,0);

f = zeros(n,1);


%Anear = zeros(d * (d-1), n);
Anear = [];
for i = 1:d
    for j = 1:d
        if i ~= j
            con_row = zeros(1,n);
            con_row(1, 2*i-1) = 1;
            con_row(1, 2*j-1) = -1;
            Anear = [Anear; con_row];
        end
    end
end
%bnear = ones(size(Anear, 1), 1);
%near_poly = Polyhedron(Anear, bnear);

%Anearv = zeros(d * d, n);
Anearv = [];
for i = 1:d
    for j = 1:d
        if i ~= j
            con_row = zeros(1,n);
            con_row(1, 2*i) = 1;
            con_row(1, 2*j) = -1;
            Anearv = [Anearv; con_row];
        end
    end
end
%bnearv = ones(size(Anearv, 1), 1);
%near_polyv = Polyhedron(Anearv, bnearv);

Apos = [-kron(eye(d), [1, 0]);
        kron(eye(d), [1, 0]);
        kron(eye(d), [0, 1]);
        -kron(eye(d), [0, 1])];
bpos = [Mh * ones(d, 1);
        Mh * ones(d, 1);
        Mv * ones(d,1);
        Mv * ones(d,1)];

%X = (kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n)) &  (delta * near_poly); % the safe state space
%X = kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n);
X = Polyhedron(Apos, bpos);% & kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n);
U = uH * Polyhedron.unitBox(1); % the input space
W = wH * Polyhedron.unitBox(d); % the disturbance space
%W = Polyhedron.fullSpace(0);
V = Polyhedron.fullSpace(0);

sys = LTISys(A,B,E,C,G);

%% j

P = [Anear;Anearv;Apos];
q0 = zeros(size(P,1),1);
%r0 = zeros(n,1);
%q0 = q1s;
%q0 = qqq;

%uset = Polyhedron.fullSpace(n) * polyhedronPower(W,T) * polyhedronPower(V,T);
%aset = polyhedronPower(U,T) * polyhedronPower(X,T+1);
uset = polyhedronProduct({Polyhedron.fullSpace(n), polyhedronPower(W,T), polyhedronPower(V,T)});
%aset = polyhedronProduct({polyhedronPower(U,T), polyhedronPower(X,T+1)});
aset = polyhedronProduct({polyhedronPower(X,T+1), polyhedronPower(U,T)});

cf = ControlFamily(sys, T, 1);

options = sdpsettings('verbose', 2, 'solver', 'gurobi'); % options for the LP solver

%c = rand(size(P,1),1)+0.05;
c = ones(size(P,1),1);

tic
[q, converged, infeasible] = insideOutRecurrenceOpt(sys,cf,uset,aset, P, q0, c, 10, 0.05, options);
toc
%[q, converged, infeasible] = insideOutRecurrence(sys,cf,uset,aset, P, q0, c, 10, 0.05, options);

converged
infeasible

%%

nsuc = 0;
nfail = 0;
suctime = 0;
failtime= 0;

%%

XF = Polyhedron('A',P,'b',q);

for i = 1:400
    [i,nsuc]

    %x0 = XF.chebyCenter.x;

    x0_test = 1*randn(2*d,1);
    tic;
    isReach = isAffineBackwardsReachable(sys,T,x0_test,X*U,XF,W,V,options);
    t = toc;
    if isReach == 1
        nsuc = nsuc + 1;
        suctime = suctime + t;
    else
        nfail = nfail + 1;
        failtime = failtime + t;
    end
end

