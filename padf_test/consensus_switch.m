%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
T = 8; % number of steps
d = 6; % number of UAVs
n = 2 * d; % dimension of state space
m = 1; % dimension of input space
l = d; % dimension of disturbance space
ns = 2; % number of switching modes

Kh = 1;
Kv = 1;

topos = cell(ns,1);

grid_topo = createGridGraph(2,3);

topo1 = grid_topo;
topo1(4,5) = 0;
topo1(5,4) = 0;
%topo1(2,3) = 0;
%topo1(3,2) = 0;


%{
topo1(1,4) = 0;
topo1(4,1) = 0;
%}
%{
topo1(1,4) = 0;
topo1(4,1) = 0;
topo1(3,6) = 0;
topo1(6,3) = 0;
topo1(4,7) = 0;
topo1(7,4) = 0;
topo1(6,9) = 0;
topo1(9,6) = 0;
%}
topos{1} = topo1;
topo2 = grid_topo;
topo2(1,2) = 0;
topo2(1,2) = 0;
topo2(5,6) = 0;
topo2(6,5) = 0;
%{
topo2(2,3) = 0;
topo1(3,2) = 0;
%}
%{
topo1(1,2) = 0;
topo1(2,1) = 0;
topo1(2,3) = 0;
topo1(3,2) = 0;
topo1(7,8) = 0;
topo1(8,7) = 0;
topo1(8,9) = 0;
topo1(9,8) = 0;
%}
topos{2} = topo2;
%topos{1} = topos{2};
%topos{1} = grid_topo;
topos{2} = grid_topo;


Mh = 8;
offset = 0;
Mv = 5;
uH = 2;
wH = 0.0005;

Mh = 20;
offset = 0;
Mv = 10;
uH = 3;
wH = 0.05;


A = cell(ns,1);
B = cell(ns,1);
E = cell(ns,1);
f = cell(ns,1);

for moded = 1:ns
    topo = topos{moded};
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

    A{moded} = Aol + Acl;
    %B{moded} = [zeros(floor(n/2)-1, 1); 0;1;zeros(ceil(n/2)-1, 1)];
    B{moded} = repmat([0;1],d,1);
    E{moded} = kron(eye(d),[0;1]);
    f{moded} = zeros(n,1);
end

%B{i} = [0;1;zeros(n-2, 1)];
%B = [zeros(floor(n/2)-1, 1); 0;1;zeros(ceil(n/2)-1, 1)];
%B{i} = repmat([0;1], d,1);
%B{i} = kron(eye(d), [0;1]);

%E = kron(eye(d),[0;1]);
%E = zeros(n,0);
%E{i} = zeros(n,0);


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

X = Polyhedron(Apos, bpos);% & kron(eye(d), [Mh, 0; 0, Mv]) * Polyhedron.unitBox(n);
U = uH * Polyhedron.unitBox(1); % the input space
W = wH * Polyhedron.unitBox(d); % the disturbance space
XU = X*U;

XUmap = cell(T,1);
WSigmamap = cell(T,ns);
for t = 1:T
        XUmap{t} = XU;
    for moded = 1:ns
        WSigmamap{t,moded} = W;
    end
end

sys = LTVSSys(T,A,B,E,f,XUmap,X,WSigmamap);
%% j

P = [Anear;
     Anearv;
     Apos];

%{
q = [ones(size(Anear,1),1);
     ones(size(Anearv,1),1);
     ones(size(Apos,1),1)];
%}
     %use q from nonswitch consensus inside out
Omega = Polyhedron('A',P,'b',q);

options = sdpsettings('verbose', 2, 'solver', 'gurobi'); % options for the LP solver

1
[isRecurrent, affineController] = testAffineRecurrence(sys,Omega,T);
isRecurrent

%%

close all
figure(1);
hold on;
xlim([1,T
]);
reps = 1;
xlabel('Time t');
ylabel('UAV Height h(t)');
title('UAV Heights under recurrent controller');

for i = 1:6
    x0_test = .2*randn(2*d,1);
    useq = [];
    xseq = [x0_test];
    x = x0_test;
    u = zeros(m,1);
            
    if Omega.contains(x0_test)
        x0 = x0_test;
        for rep = 1:reps
            if rep > 1
                x0 = xseq(:,end);
            end
            w_seq = [];
            modes = [];
            seq = '';
            for t = 1:T
                wind = randi(size(W.V,1));
                mode = randi(ns);
                w_mult = rand()/2;
                w_seq = [w_seq; w_mult*W.V(wind,:)'];
                modes = [modes, mode];
            end

            for t = 1:(T-1)
                seq = LTVSSys.getSequenceFromModes(modes(1:t));
                u = affineController.Kx_map(seq) * x0 + affineController.Kw_map(seq)*w_seq+affineController.uc_map(seq);
                x = A{modes(t)} * x + B{modes(t)} * u + E{modes(t)} * w_seq((t-1)*l+(1:l)) + f{modes(t)};
                useq = [useq, u];
                xseq = [xseq, x];
            end
        end
        trange = 1:(1+(reps*(T-1)));
        colors = {'r','g','b',[0.9,0.9,0],[0.9,0,0.6],'k'};
        for uav = 1:d
            plot(trange,xseq(1+(uav-1)*2,:),'color',colors{uav});
        end
    end    
end

%%

[constraints, x0,lambda,x0_c,Kw_maps,uc_maps] = computeAffineBackwardsReachableSet(sys,Omega);

nsuc = 0;
nfail = 0;
suctime = 0;
failtime= 0;
for i = 1:400
    [i,nsuc]
    
    x0_test = 4*randn(2*d,1);
    tic;
    diagnostics = optimize([constraints; x0==x0_test], [], options);
    t = toc;
    if diagnostics.problem == 0
        nsuc = nsuc + 1;
        suctime = suctime + t;
    else
        nfail = nfail + 1;
        failtime = failtime + t;
    end
end
