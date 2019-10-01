function [diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(Omega, pslsys, N, options)
% evalSwitchedInvariance determines if the given switched system and seed
%   polytope generate a controlled invariant set robust to switching
%   
%   More precisely this function uses a linear program to determine if
%   there is a control policy that drives all points of Omega into Omega
%   in exactly N steps. This policy is for a fixed switching history,
%   linear in the initial state and disturbance history. This guarantees
%   that conv(\bigcup_{i=1}^N Pre^i(Omega)) is an invariant set.
%
%   Omega - seed set polytope
%   pslsys - the polytopic switching linear system
%   N - the number of steps to consider
%   options - options for the Yalmip LP solver
%
%   diagnostics - results from the Yalmip solver
%   sequences - all possible switching sequences
%   Kx_map - the map from switching histories to initial state feedback
%   Kw_map - the map from switching histories to disturbance feedback
%   uc_map - the map from switching histories to input offset
%

% compute all possible switching sequences and the corresponding system quantities
[sequences, A_seq, f_seq] = computeSwitchedSequence(pslsys.A, pslsys.f, N);

n = size(pslsys.A{1},2); % state dim
m = size(pslsys.B{1},2); % input dim
l = size(pslsys.E{1},2); % disturbance dim

nUcon = size(pslsys.U.A,1); % number of constraints defining U
nWcon = size(pslsys.W.A,1); % number of constraints defining W

% Each partial switching sequence corresponds to a target constraint either
% enforcing safety or reaching the target of Omega
% These constraints are of the form
%       Gx * x + Gu * u + Gw * w <= g
%
% The following variables define maps from partial sequences to the
% constraint matrices/vectors. They are initialized with the constraint
% corresponding to the empty string ''.
target_G_map_x = containers.Map({''},{pslsys.X.A});
target_G_map_u = containers.Map({''},{zeros(size(pslsys.X.A,1),N * m)});
target_G_map_w = containers.Map({''},{zeros(size(pslsys.X.A,1),N * l)});
target_g_map = containers.Map({''},{pslsys.X.b});

% We define these constraints for all lengths of sequences from 1 to N
for i = 1:N
    % For all but the last step, the target set is simply the safe set
    if i == N
        target = Omega;
    else
        target = pslsys.X;
    end
    len_sequences = sequences{i+1};
    for j = 1:size(len_sequences, 1)
        seq = len_sequences{j};
        [Gx, Gu, Gw, g] = computeSwitchedTargetConstraint(pslsys, A_seq, f_seq, target, seq, N);
        target_G_map_x(seq) = Gx;
        target_G_map_u(seq) = Gu;
        target_G_map_w(seq) = Gw;
        target_g_map(seq) = g;
    end
end

% clear LP solver variables
yalmip('clear')
constraints = [];

% These quantities define the polytop Omega x W^N = {z | Hbar * z <= hbar}
Hbar = blkdiag(Omega.A,kron(eye(N),pslsys.W.A));
hbar = [Omega.b;repmat(pslsys.W.b,N,1)];

% We will search for control policies that are affine in the initial state
% and the disturbance history given the switching history.
% These containers map partial sequences to the appropriate decision
% variable corresponding to the affine control law
%       u = Kx * x + Kw * w + uc
Kx_map = containers.Map('KeyType','char','ValueType','any');
Kw_map = containers.Map('KeyType','char','ValueType','any');
uc_map = containers.Map('KeyType','char','ValueType','any');
for i = 1:N
    len_sequences = sequences{i};
    for j = 1:size(len_sequences, 1)
        seq = len_sequences{j};
        Kx = sdpvar(m, n, 'full');
        Kx_map(seq) = Kx;
        % Kw maps the known disturbance history to an input, but is padded
        % with zeros to be a function of the entire history
        Kw = [sdpvar(m, (i-1) * l, 'full'), zeros(m, (N - i + 1) * l)];
        Kw_map(seq) = Kw;
        ux = sdpvar(m,1,'full');
        uc_map(seq) = ux;
    end
end

% For each switching sequence of full length we add appropriate constraints
total_sequences = sequences{N+1};
for j = 1:size(total_sequences, 1)

    sequence = total_sequences{j};
    
    % We build the constraint matrices in block form
    target_G_x = cell(N+1,1);
    target_G_u = cell(N+1,1);
    target_G_w = cell(N+1,1);
    target_g = cell(N+1,1);
    
    split_seq = split(sequence,',');
    partial_seq = '';
    
    % The blocks of the control decision variable matrices are already
    % allocated in the container maps Kx_map, Kw_map, uc_map, but in order
    % to build the matrix for each sequence from these blocks, we must
    % reallocate variables even though they will be overwritten.
    Kx_inst = sdpvar(N * m, n, 'full');
    Kw_inst = sdpvar(N * m, N * l, 'full');
    uc_inst = sdpvar(N * m, 1, 'full');
    for i = 1:(N+1)
        target_G_x{i} = target_G_map_x(partial_seq);
        target_G_u{i} = target_G_map_u(partial_seq);
        target_G_w{i} = target_G_map_w(partial_seq);
        target_g{i} = target_g_map(partial_seq);
        
        if i < N+1
            Kx_inst(m * (i-1) + (1:m), :) = Kx_map(partial_seq);
            Kw_inst(m * (i-1) + (1:m), :) = Kw_map(partial_seq);
            uc_inst(m * (i-1) + (1:m), :) = uc_map(partial_seq);
            if i == 1
                partial_seq = split_seq{i};
            else
                partial_seq = strcat(partial_seq, ',', split_seq{i});
            end
        end
    end
    
    % combine the blocks into constraint matrices with the correct
    % dimension by padding as necessary
    Gx = [cell2mat(target_G_x); zeros(N * (nUcon + nWcon), n)];
    Gu = [cell2mat(target_G_u); kron(eye(N),pslsys.U.A); zeros(N * nWcon, N * m)];
    Gw = [cell2mat(target_G_w); zeros(N * nUcon, N * l); kron(eye(N),pslsys.W.A)];
    g = [cell2mat(target_g);  repmat(pslsys.U.b, N, 1); repmat(pslsys.W.b, N, 1)];
    
    gu_offset = Gu * uc_inst;
    
    T = sdpvar(size(Gx,1),  size(Hbar,1), 'full');   
    
    % we add the polytope inclusion constraints for this sequence to our
    % existing list
    constraints = [constraints;
        T*Hbar == [Gx + Gu * Kx_inst, Gu * Kw_inst + Gw];
        T*hbar <= g - gu_offset;
        T >= 0;
    ];

end

% finally we solve the linear feasibility program with no objective
diagnostics = optimize(constraints, [], options);

% the control policy determined can be accessed through the maps
% Kx_map, Kw_map, uc_map indexed by elements of sequences

end