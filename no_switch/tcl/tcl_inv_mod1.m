
Cf = eye(Kp);

Ef = Bf;
Gf = zeros(Kp,0);

W = 0.5 * Polyhedron.unitBox(m);
V = Polyhedron.fullSpace(0);

sys = LTISys(Af,Bf,Ef,Cf,Gf);

% The horizon considered
Tpre = 0;
Tcyc = 8;
T = Tpre + Tcyc;

epsilon = 0.01;

% This polyhedron models that constraints for the upper and lower bound on
% the number of active loads, the nonnegativity requirements, and the total
% number of loads.
off_con = zeros(1,Kp);
off_con(off_ind) = 1;
on_con = zeros(1,Kp);
on_con(on_ind) = 1;
n_off = length(off_ind);
n_on = length(on_ind);
powerConstraintA = ...
    [off_con;
    -eye(n_off), zeros(n_off,Kp-n_off);
    on_con;
    zeros(n_on,Kp-n_on), -eye(n_on);
    ];
powerConstraintb = ...
    [N - lowerBound;
    zeros(n_off,1);
    upperBound;
    zeros(n_on,1);
    ];

% This polyhedron models the constraints that the number of loads to switch
% from m1 to m2 cannot be more than the number of loads at m1, and that the
% number of loads switched must be nonnegative.
input_x_con = eye(Kp);
input_x_con = input_x_con(input_ind,:);
jointInputConstraintA_x = [
    -input_x_con;
    zeros(m,Kp)
    ];
jointInputConstraintA_u = [
    eye(m);
    -eye(m)
    ];
jointInputConstraintb = [
    zeros(m,1);
    zeros(m,1)
    ];

jointInputConstraintAHor = ...
    [kron(eye(T), jointInputConstraintA_x), ...
     zeros(T*size(jointInputConstraintA_x,1),Kp), ...
     kron(eye(T), jointInputConstraintA_u)
     ];
jointInputConstraintbHor = ...
    repmat(jointInputConstraintb, T, 1);
 
%% Invariant set computations
%permMat = [blkdiag(kron(eye(T), [eye(2*K);zeros(2*K)]), eye(2*K)), ...
%           [kron(eye(T), [zeros(2*K); eye(2*K)]); zeros(2*K,2*K*T)]];
1.5
A_ad = [[kron(eye(T+1),powerConstraintA), zeros((T+1)*size(powerConstraintA,1), m*T)];
        jointInputConstraintAHor;
        %cycConstraintA
        ];
b_ad = [repmat(powerConstraintb, T+1,1);
        jointInputConstraintbHor;
        %cycConstraintb
        ];
3
cf = ControlFamily(sys, T, 1, 0);
4
options = sdpsettings('verbose', 0, 'solver', 'gurobi'); % options for the LP solver

%%

%x0 = N/Kp*ones(Kp,1);
%x0 = N/12*[1;2;2;2;1;1;1;1;1];

%x0 = N/12*[1;1;2;2;1;1;1;1;1];
x0 = ones(Kp,1);
x0(9) = 2;
x0(8) = 2;
x0(6) = 2;
x0 = x0 * N/sum(x0);
%Omega = x0 + 100*W;
%Omega = Polyhedron('V',x0');
%Omega = Polyhedron('A', powerConstraintA, 'b', powerConstraintb);

Omega = 20*Polyhedron.unitBox(Kp);
Omega = Polyhedron('A',Omega.A,'b',Omega.b,'Ae',ones(1,Kp),'be',0);
Omega = x0+Omega;
A_Omega = [Omega.A; Omega.Ae; -Omega.Ae];
b_Omega = [Omega.b; Omega.be; -Omega.be];

A_ad = [A_ad; 
        zeros(size(A_Omega,1),Kp*T),A_Omega,zeros(size(A_Omega,1),m*T)];
b_ad = [b_ad; b_Omega];
%b_ad = b_ad + epsilon *ones(size(b_ad,1),1);

cf = ControlFamily(sys, T, 1,1);
options = sdpsettings('verbose', 2, 'solver', 'gurobi'); % options for the LP solver

uset = polyhedronProduct({Omega, polyhedronPower(W,T), polyhedronPower(V,T)});
aset = Polyhedron('A',A_ad,'b',b_ad);

tic
diagnostics = computeController(sys,cf,uset,aset, options);
toc

diagnostics.problem

%%

options = sdpsettings('verbose', 0, 'solver', 'gurobi'); % options for the LP solver
T=T;
nsuc = 0;
nfail = 0;
suctime = 0;
failtime= 0;
XU = Polyhedron('A',[powerConstraintA, zeros(size(powerConstraintA, 1),m); jointInputConstraintA_x,jointInputConstraintA_u], 'b', [powerConstraintb; jointInputConstraintb]);
XF = Omega;
for i = 1:100
    [i,nsuc]
    xtest = x0 + 100*randn(Kp,1);
    xtest = xtest / sum(xtest) * N;
    tic;
    isReach = isAffineBackwardsReachable(sys,T,xtest,XU,XF,W,V,options);
    t = toc;
    if isReach
        nsuc = nsuc + 1;
        suctime = suctime + t;
    else
        nfail = nfail + 1;
        failtime = failtime + t;
    end
end

%{
x0 = N/2/K*ones(2*K,1);
5
diag = computeControllerForPointNoDist(sys,cf,x0,aset, options);
%}
