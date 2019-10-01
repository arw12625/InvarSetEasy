%% system parameters

% The switched affine dynamics are defined by
%       x+ = A{sigma} * x + B{sigma} * u + E{sigma} * w + f{sigma}
N = 5; % number of steps
n = 2; % dimension of state space
m = 2; % dimension of input space
l = 1; % dimension of disturbance space
ns = 2; % number of switching modes

A = cell(ns,1); % A matrices for each mode
A{1} = 1.4 * eye(n);
A{2} = 1.2 * eye(n);

B = cell(ns,1); % B matrices for each mode
B{1} = [1,0;0.1,1];
B{2} = [1,1;0,1];

E = cell(ns,2); % E matrices for each mode
E{1} = eye(n);
E{2} = eye(n);

f = cell(ns,1); % f vectors for each mode
f{1} = zeros(n,1);
f{2} = zeros(n,1);

Omega = 1.7 * Polyhedron.unitBox(n); % the seed set
X = 10 * Polyhedron.unitBox(n); % the safe state space
U = 1 * Polyhedron.unitBox(n); % the input space
W = 0.1 * Polyhedron.unitBox(2); % the disturbance space

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f); % class representing the system

yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver

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