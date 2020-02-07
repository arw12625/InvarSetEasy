function [diagnostics, affineController] = computeAffineController(sys, initialConditions, admissibleTrajectories, yalmipOptions)
    % computeAffineController constructs an affine controller if possible
    % that drives all of the initial conditions into admissible
    % trajectories through the system dynamics.
    % This is an affine disturbance feedback controller.
    %
    %   sys - an LTVSSys representing the system
    %   initialConditions - initial conditions (see computeInitialConditions)
    %   admissibleTrajectories - admissible lifted trajectories (see admissibleTrajectories)
    %   yalmipOptions - options for the linear program solver
    %
    %   diagnostics - output from the linear program solver
    %   affineController - a structure encompassing the affine controller
    %                   if found. The controller ac maps a switching sequence
    %                   sigma, initial state x0, and disturbance sequence w
    %                   to the input 
    %                   u = ac.Kx_map(sigma) * x0 + ac.Kw_map(sigma) * w + ac.uc_map(sigma)
    %
assert(initialConditions.horizon == admissibleTrajectories.horizon);
horizon = initialConditions.horizon;

n = sys.n;
m = sys.m;
l = sys.l;

% clear LP solver variables
yalmip('clear')
constraints = [];

% We will search for control policies that are affine in the initial state
% and the disturbance history given the switching history.
% These containers map partial sequences to the appropriate decision
% variable corresponding to the affine control law
%       u = Kx * x + Kw * w + uc
Kx_map = containers.Map('KeyType','char','ValueType','any');
Kw_map = containers.Map('KeyType','char','ValueType','any');
uc_map = containers.Map('KeyType','char','ValueType','any');
for i = 1:horizon
    len_sequences = sys.sequences{1,i};
    for j = 1:size(len_sequences, 1)
        seq = len_sequences{j};
        Kx = sdpvar(m, n, 'full');
        Kx_map(seq) = Kx;
        % Kw maps the known disturbance history to an input, but is padded
        % with zeros to be a function of the entire history
        Kw = [sdpvar(m, (i-1) * l, 'full'), zeros(m, (horizon - i + 1) * l)];
        Kw_map(seq) = Kw;
        ux = sdpvar(m,1,'full');
        uc_map(seq) = ux;
    end
end


% For each switching sequence of full length we add appropriate constraints
total_sequences = sys.sequences{1,horizon+1};
for j = 1:size(total_sequences, 1)

    sequence = total_sequences{j};
    modes = LTVSSys.getModesFromSequence(sequence);
    
    % The blocks of the control decision variable matrices are already
    % allocated in the container maps Kx_map, Kw_map, uc_map, but in order
    % to build the matrix for each sequence from these blocks, we must
    % reallocate variables even though they will be overwritten.
    Kx_inst = sdpvar(horizon * m, n, 'full');
    Kw_inst = sdpvar(horizon * m, horizon * l, 'full');
    uc_inst = sdpvar(horizon * m, 1, 'full');
    for t = 1:horizon
        prefix = LTVSSys.getSequenceFromModes(modes(1:(t-1)));
        Kx_inst(m * (t-1) + (1:m), :) = Kx_map(prefix);
        Kw_inst(m * (t-1) + (1:m), :) = Kw_map(prefix);
        uc_inst(m * (t-1) + (1:m), :) = uc_map(prefix);
    end
    
    Gx = admissibleTrajectories.AxMap(sequence);
    Gu = admissibleTrajectories.AuMap(sequence);
    Gw = admissibleTrajectories.AwMap(sequence);
    g = admissibleTrajectories.bMap(sequence);
    
    % These quantities define the polytop Omega x W^N = {z | Hbar * z <= hbar}
    initPoly = initialConditions.initPolyMap(sequence);
    Hbar = [initPoly.A; initPoly.Ae; -initPoly.Ae];
    hbar = [initPoly.b; initPoly.be; -initPoly.be];
    
    T = sdpvar(size(Gx,1),  size(Hbar,1), 'full');   
    
    % we add the polytope inclusion constraints for this sequence to our
    % existing list
    constraints = [constraints;
        T*Hbar == [Gx + Gu * Kx_inst, Gu * Kw_inst + Gw];
        T*hbar <= g - Gu * uc_inst;
        T >= 0;
    ];
end

% finally we solve the linear feasibility program with no objective
diagnostics = optimize(constraints, [], yalmipOptions);

%(Polyhedron('H',[Hbar, hbar]) < Polyhedron('H',[value([Gx + Gu * Kx_inst, Gu * Kw_inst + Gw]), value(g - gu_offset)]))
%plot(Polyhedron(value([Gx + Gu * Kx_inst, Gu * Kw_inst + Gw]), value(g - gu_offset)))

affineController.Kx_map = containers.Map('KeyType','char','ValueType','any');
affineController.Kw_map = containers.Map('KeyType','char','ValueType','any');
affineController.uc_map = containers.Map('KeyType','char','ValueType','any');
for i = 1:horizon
    len_sequences = sys.sequences{1,i};
    for j = 1:size(len_sequences, 1)
        seq = len_sequences{j};
        affineController.Kx_map(seq) = value(Kx_map(seq));
        affineController.Kw_map(seq) = value(Kw_map(seq));
        affineController.uc_map(seq) = value(uc_map(seq));
    end
end
affineController.horizon = horizon;

end