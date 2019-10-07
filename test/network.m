%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 6; % number of steps
n = 4; % dimension of state space
m = 1; % dimension of input space
l = 4; % dimension of disturbance space
ns = 2; % number of switching modes

A = cell(ns,1); % A matrices for each mode
K = 0.5;
A{1} = [1-K, K/2, 0, K/2;
        K/2, 1-K, K/2, 0;
        0, K/2, 1-K, K/2;
        K/2, 0, K/2, 1-K];
A{1} = [1-K, 0, 0, K;
        0, 1-K, K, 0;
        0, K/2, 1-K, K/2;
        K/2, 0, K/2, 1-K];
A{2} = [1-K, K, 0, 0;
        K/2, 1-K, K/2, 0;
        0, K/2, 1-K, K/2;
        0, 0, K, 1-K];

B = cell(ns,1); % B matrices for each mode
B{1} = [1;1;1;1];
B{2} = [1;1;1;1];

E = cell(ns,2); % E matrices for each mode
E{1} = eye(n);
E{2} = eye(n);

f = cell(ns,1); % f vectors for each mode
f{1} = zeros(n,1);
f{2} = zeros(n,1);

Xm = 10;
Um = 10;
%Wm = 0.125;
Wm = 0.15;
delta = 2;
Acon = [1, -1, 0, 0;
        1, 0, -1, 0;
        1, 0, 0, -1;
        0, 1, -1, 0;
        0, 1, 0, -1;
        0, 0, 1, -1;];
Acon = [Acon; -Acon];
bcon = delta * ones(size(Acon, 1), 1);
near_poly = Polyhedron(Acon, bcon);

X = (Xm * Polyhedron.unitBox(n)) & near_poly; % the safe state space
U = Um * Polyhedron.unitBox(m); % the input space
W = (Wm * Polyhedron.unitBox(l)); % the disturbance space

Omega = (0.5 * Xm * Polyhedron.unitBox(n)) & (0.75 * near_poly); % the seed set

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f); % class representing the system

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver

%%

% determine if the seed set generates an invariant set with respect to the
% system dynamics
[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(Omega, pslsys, N, yalmipOptions);


if diagnostics.problem == 0
    disp('Omega generates an invariant set');
    %{
    diagnostics.problem
    for i = 1:(size(sequences, 1) - 1)
        len_seq = sequences{i};
        for j = 1:size(len_seq, 1)
            seq = len_seq{j}
            value(Kx_map(seq))
            value(Kw_map(seq))
        end
    end
    %}
end