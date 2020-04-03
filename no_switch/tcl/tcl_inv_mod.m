
C = eye(2 * K);

% No disturbance
E = zeros(2 * K, 0);
G = zeros(2 * K, 0); 

W = Polyhedron.fullSpace(0);
V = Polyhedron.fullSpace(0);

% Disturbance

E = B;
G = zeros(2*K,0);

W = 0.5 * Polyhedron.unitBox(2*(K-1));
V = Polyhedron.fullSpace(0);

sys = LTISys(A,B,E,C,G);

% The horizon considered
Tpre = 0;
Tcyc = 6;
T = Tpre + Tcyc;

epsilon = .1;

% This polyhedron models that constraints for the upper and lower bound on
% the number of active loads, the nonnegativity requirements, and the total
% number of loads.
powerConstraintA = ...
    [ones(1,K), zeros(1,K);
    -eye(K), zeros(K);
    zeros(1,K), ones(1,K);
    zeros(K), -eye(K);
    ];
powerConstraintb = ...
    [N - lowerBound;
    zeros(K,1);
    upperBound;
    zeros(K,1);
    ];

% This polyhedron models the constraints that the number of loads to switch
% from m1 to m2 cannot be more than the number of loads at m1, and that the
% number of loads switched must be nonnegative.
jointInputConstraintA_x = [-eye(K-1), zeros(K-1,K+1);
     zeros(K-1,2*K);
     zeros(K-1,K+1), -eye(K-1);
     zeros(K-1,2*K)];
jointInputConstraintA_u = [eye(K-1), zeros(K-1);
      -eye(K-1), zeros(K-1);
      zeros(K-1), eye(K-1);
      zeros(K-1), -eye(K-1)];
jointInputConstraintb = [zeros(K-1,1);
    zeros(K-1,1);
    zeros(K-1,1);
    zeros(K-1,1)];
%jointInputConstraintb = epsilon*ones(4*K,1);

jointInputConstraintAHor = ...
    [kron(eye(T), jointInputConstraintA_x), ...
     zeros(T*size(jointInputConstraintA_x,1),2*K), ...
     kron(eye(T), jointInputConstraintA_u)
     ];
jointInputConstraintbHor = ...
    repmat(jointInputConstraintb, T, 1);
 
%{
cycConstraintA = ...
    [zeros(2*K, 2*K*Tpre), eye(2*K), zeros(2*K, 2*K*(Tcyc-1)), -eye(2*K), zeros(2*K,2*K*(Tpre+Tcyc));
    zeros(2*K, 2*K*Tpre), eye(2*K), zeros(2*K, 2*K*(Tcyc-1)), -eye(2*K), zeros(2*K,2*K*(Tpre+Tcyc))];
cycConstraintb = zeros(2*2*K,1);
%}


%% Invariant set computations
%permMat = [blkdiag(kron(eye(T), [eye(2*K);zeros(2*K)]), eye(2*K)), ...
%           [kron(eye(T), [zeros(2*K); eye(2*K)]); zeros(2*K,2*K*T)]];
1.5
A_ad = [[kron(eye(T+1),powerConstraintA), zeros((T+1)*size(powerConstraintA,1), 2*(K-1)*T)];
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

x0 = N/2/K*ones(2*K,1);
%Omega = x0 + 100*W;
%Omega = Polyhedron('V',x0');
%Omega = Polyhedron('A', powerConstraintA, 'b', powerConstraintb);

Omega = 10*Polyhedron.unitBox(2*K);
Omega = Polyhedron('A',Omega.A,'b',Omega.b,'Ae',ones(1,2*K),'be',0);
A_Omega = [Omega.A; Omega.Ae; -Omega.Ae];
b_Omega = [Omega.b; Omega.be; -Omega.be];
Omega = x0+Omega;

A_ad = [A_ad; 
        zeros(size(A_Omega,1),2*K*T),A_Omega,zeros(size(A_Omega,1),2*(K-1)*T)];
b_ad = [b_ad; b_Omega+epsilon];
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

%{
XU = Polyhedron('A',[jointInputConstraintA_x,jointInputConstraintA_u], 'b', jointInputConstraintb);
XF = Omega;
isReach = isAffineBackwardsReachable(sys,T,x0,XU,XF,W,V,options)
%}

%{
x0 = N/2/K*ones(2*K,1);
5
diag = computeControllerForPointNoDist(sys,cf,x0,aset, options);
%}
