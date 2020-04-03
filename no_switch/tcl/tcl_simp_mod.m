%% System definitions

% This script investigates invariant set computations for a TCL system


% Number of discrete states
K = 5;

% These adjacency matrices model the simple case where in one mode the
% discrete temperature increments by one while in the other mode it 
% decrements by one at each time step.
N0_adjacency = zeros(K,K);
N0_adjacency(2:end,1:end-1) = eye(K-1);
N0_adjacency(end,end) = 1;

N1_adjacency = zeros(K,K);
N1_adjacency(1:end-1,2:end) = eye(K-1);
N1_adjacency(1,1) = 1;

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
A = blkdiag(N0_adjacency, N1_adjacency);
B = [-N0_adjacency, N0_adjacency; N1_adjacency, -N1_adjacency];

Kf = diag([zeros(K-2,1);1;0;0;1;zeros(K-2,1)]);
Af = A+B*Kf;

state_ind = find([0;0;ones(K-3,1);0;0;ones(K-3,1);0;0]);
input_ind = find([0;0;ones(K-4,1);0;0;0;0;ones(K-4,1);0;0]);

Amod = Af(state_ind,state_ind);
Bmod = B(state_ind,input_ind);
A = Amod;
B = Bmod;
K = K - 3;

% System constraints as polyhedra
N = 10000; % number of loads
upperBound = 7000; % upper bound on the number of loads that are on
lowerBound = 3000; % lower bound on the number of loads that are on
%upperBound = N; % upper bound on the number of loads that are on
%lowerBound = 0; % lower bound on the number of loads that are on