function [poly, x0, constraints] = computeLiftedAffineSat(sys, admissibleTrajectories)
    % computeLiftedAffinePre computes a Polyhedron consisting of pairs of
    % an initial state and an affine controller that drives that state
    % along an admissible trajectory
    
    horizon = admissibleTrajectories.horizon;
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
    %       u = Kw * w + uc
    Kw_map = containers.Map('KeyType','char','ValueType','any');
    uc_map = containers.Map('KeyType','char','ValueType','any');
    for i = 1:horizon
        len_sequences = sys.sequences{1,i};
        for j = 1:size(len_sequences, 1)
            seq = len_sequences{j};
            % Kw maps the known disturbance history to an input, but is padded
            % with zeros to be a function of the entire history
            Kw = [sdpvar(m, (i-1) * l, 'full'), zeros(m, (horizon - i + 1) * l)];
            Kw_map(seq) = Kw;
            ux = sdpvar(m,1,'full');
            uc_map(seq) = ux;
        end
    end

    yalmipOptions = sdpsettings('verbose', 1); % options for the LP solver
    total_sequences = sys.sequences{1,horizon+1};
    constraints = [];
    x0 = sdpvar(n,1,'full');
    
    
    initialConditions = computeInitialConditions(sys, Polyhedron.fullSpace(n), horizon, 0, 1);
    for j = 1:size(total_sequences, 1)

        sequence = total_sequences{j};
        modes = LTVSSys.getModesFromSequence(sequence);

        % The blocks of the control decision variable matrices are already
        % allocated in the container maps Kx_map, Kw_map, uc_map, but in order
        % to build the matrix for each sequence from these blocks, we must
        % reallocate variables even though they will be overwritten.
        Kw_inst = sdpvar(horizon * m, horizon * l, 'full');
        uc_inst = sdpvar(horizon * m, 1, 'full');
        for t = 1:horizon
            prefix = LTVSSys.getSequenceFromModes(modes(1:(t-1)));
            Kw_inst(m * (t-1) + (1:m), :) = Kw_map(prefix);
            uc_inst(m * (t-1) + (1:m), :) = uc_map(prefix);
        end

        Gx = admissibleTrajectories.Gxmap(sequence);
        Gu = admissibleTrajectories.Gumap(sequence);
        Gw = admissibleTrajectories.Gwmap(sequence);
        g = admissibleTrajectories.gmap(sequence);

        % These quantities define the polytop Omega x W^N = {z | Hbar * z <= hbar}
        Hbar = initialConditions.Hmap(sequence);
        hbar = initialConditions.hmap(sequence);


        T = sdpvar(size(Gx,1),  size(Hbar,1), 'full');   

        %size(x0)
        %size(Kw_inst)
        %size(uc_inst)
        %size(T)
        
        % we add the polytope inclusion constraints for this sequence to our
        % existing list
        constraints = [constraints;
            T*Hbar == Gu * Kw_inst + Gw;
            T*hbar <= g - Gx * x0 - Gu * uc_inst;
            T >= 0;
        ];
    end
    
    poly = Polyhedron(constraints);
    
end


