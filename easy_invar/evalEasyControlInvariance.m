function [isInvar,diagnostics,Kx,Kw,u0] = evalEasyControlInvariance(Omega, plsys, N, useDistFeedback, useAffineControl, yalmipOptions)
%evalEasyControlInvariance Determine if Omega \subset pre^N(Omega) in which
%case conv(\cup_{i=1}^N pre^N(Omega)) is a controlled invariant set
%   More specifically, the program searches for a feedback controller that
%   drives Omega into itself in N steps. This method is sound but not
%   complete
%
%   Omega - seed set
%   plsys - polytopic linear system
%   N - number of backwards steps
%   useDistFeedback - whether or not the controller should use disturbance
%       feedback (if not then the control depends only on the initial
%       state, which yields very conservative results
%   useAffineControl - whether or not the controller has an affine term
%   yalmipOptions - options for the solver (note if the system has
%       disturbance, the solver specified must handle bilinear constraints)
%
%   isInvar - whether Omega generates a controlled invariant set
%   diagnostics - results from the solver
%   Kx,Kw,u0 - the controller found by the algorithm
%       u = Kx x_0 + Kw w + u0
%       where x_0 is the initial state and w is the disturbance sequence
%       note that u,w are in reverse order, the last input is first

yalmip('clear')

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


[Gbarx,Gbaru,Gbarw,gbarh, gbarg] = computeLiftedPre(plsys,Omega,N,1);



Hbar = blkdiag(H,kron(eye(N),E));
hbar = [h;repmat(e,N,1)];


%Decision variables for LP
T = sdpvar(size(Gbarx,1),  size(Hbar,1), 'full');
Kx = sdpvar(N * m, n, 'full');

Kw = sdpvar(N*m, N*l,'full');
for i = 1:N
    Kw((i-1)*m + (1:(N-i+1)*m),(i-1)*l+(1:l)) = zeros((N-i+1)*m,l);
end

if ~useDistFeedback
    Kw = Kw * 0;
end

u0 = sdpvar(N*m, 1,'full');
if ~useAffineControl
    u0 = u0 * 0;
end

%{
targetInputCon = cell(1,N);
targetDistCon = cell(1,N);
stateStateCon = cell(N+1,1);
stateInputCon = cell(N+1,N);
stateDistCon = cell(N+1,N);

gbar = cell(1+(N+1)+N+N,1);

gbar{1,1} = h;
gbar{1+N+1,1} = f;
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
    gbar{1+i,1} = f;
    gbar{1+(N+1)+i,1} = g;
    gbar{1+(N+1)+N+i,1} = e;
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

Gbarx = cell2mat(Gbarx);
Gbaru = cell2mat(Gbaru);
Gbarw = cell2mat(Gbarw);

gbar = cell2mat(gbar);
%}

constraints = [
    T*Hbar == [Gbarx + Gbaru * Kx, Gbaru * Kw + Gbarw];
    T*hbar <= gbarg + gbarh - Gbaru * u0;
    T >= 0;
];

diagnostics = optimize(constraints, [], yalmipOptions);

if diagnostics.problem == 0
    isInvar = 1;
    Kx = value(Kx);
    Kw = value(Kw);
    u0 = value(u0);
else
    isInvar = 0;
    Kx = zeros(size(Kx));
    Kw = zeros(size(Kw));
    u0 = zeros(N * m,1);
end

%{
tmp = value(T*hbar - gbar);
tmp(size(H,1)+N*size(F,1) + (1:size(F,1)))
value(Kw)

check(constraints)
%}
check(constraints)
end

