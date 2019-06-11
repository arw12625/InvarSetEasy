function [Gbarx, Gbaru, Gbarw, gbarh, gbarg] = computeLiftedPre(plsys, Omega, N, incDistConst)
%computeLiftedPre Compute a polytopic representation of the lifted space
%corresponding to the N step pre of Omega
%   The computed set consists of all triples (x_0,u,w) of initial state,
%   input sequence, and disturbance sequence such that the corresponding
%   trajectory satisfies all safety,target,input,and disturbance
%   constraints over N steps.
%
%   plsys - polytopic linear system
%   Omega - target set
%   N - number of backwards steps
%   incDistConst - whether or not to include the disturbance constraints in
%       the resulting polytope
%
%   The resulting polytope is 
%       { (x_0,u,w) | Gbarx x_0 + Gbaru u + Gbarw w <= gbarh + gbarg }
%   where gbarh corresponds to the target constraints and gbarg to all
%   other constraints.
%   Note that the input and disturbance sequences u,w are specified in
%   reverse order

%Linear inequality representation of polytopes X,U,Omega
F = plsys.X.A;
f = plsys.X.b;
G = plsys.U.A;
g = plsys.U.b;
H = Omega.A;
h = Omega.b;
E = plsys.W.A;
e = plsys.W.b;

%Dimensions of system and constraints
n = plsys.n;
m = plsys.m;
l = plsys.l;

A = plsys.A;
B = plsys.B;
C = plsys.E;
%f = plsys.f;

targetInputCon = cell(1,N);
targetDistCon = cell(1,N);
stateStateCon = cell(N+1,1);
stateInputCon = cell(N+1,N);
stateDistCon = cell(N+1,N);

gbarg = cell(1+(N+1)+N+N,1);

gbarg{1,1} = zeros(size(h,1),1);
gbarg{1+N+1,1} = f;
stateStateCon{N+1,1} = F;
for i = 1:N
    targetInputCon{1,i} = H * A^(i-1) * B;
    targetDistCon{1,i} = H * A^(i-1) * C;
    stateStateCon{i,1} = F * A^(N-i+1);
    for j = i:N
        stateInputCon{j-i+1,j} = F * A^(i-1) * B;
        stateDistCon{j-i+1,j} = F * A^(i-1) * C;
        if i > 1
            stateInputCon{j,j-i+1} = zeros(size(F,1),m);
            stateDistCon{j,j-i+1} = zeros(size(F,1),l);
        end
        stateInputCon{N+1,i} = zeros(size(F,1),m);
        stateDistCon{N+1,i} = zeros(size(F,1),m);
    end
    gbarg{1+i,1} = f;
    gbarg{1+(N+1)+i,1} = g;
    gbarg{1+(N+1)+N+i,1} = e;
end

Gbarx = cell(4,1);
Gbaru = cell(4,1);
Gbarw = cell(4,1);

Gbarx{1,1} = H * A^N;
Gbaru{1,1} = cell2mat(targetInputCon);
Gbarw{1,1} = cell2mat(targetDistCon);

Gbarx{2,1} = cell2mat(stateStateCon);
Gbaru{2,1} = cell2mat(stateInputCon);
Gbarw{2,1} = cell2mat(stateDistCon);

Gbarx{3,1} = zeros(size(G,1) * N,n);
Gbaru{3,1} = kron(eye(N), G);
Gbarw{3,1} = zeros(size(G,1) * N, l * N);

Gbarx{4,1} = zeros(size(E,1) * N,n);
Gbaru{4,1} = zeros(size(E,1) * N, m * N);
Gbarw{4,1} = kron(eye(N), E);

if ~incDistConst
    Gbarx = Gbarx(1:3);
    Gbaru = Gbaru(1:3);
    Gbarw = Gbarw(1:3);
    gbarg = gbarg(1:(1+(N+1)+N));
end

Gbarx = cell2mat(Gbarx);
Gbaru = cell2mat(Gbaru);
Gbarw = cell2mat(Gbarw);

gbarg = cell2mat(gbarg);
gbarh = zeros(size(gbarg,1),1);

gbarh(1:size(h,1)) = h;

end

