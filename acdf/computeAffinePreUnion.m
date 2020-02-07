function [constraints, x0,lambda,x0_c,Kw_maps,uc_maps] = computeAffinePreUnion(base_sys, target)
% computePreUnion computes the convex hull of the union of the sets produced by iteratively
% applying the pre operator to the target set for the given number of steps
%
%
% sys - the LTVSSys representing system dynamics
% target - the polyhedron to apply pre to
% steps - the number of iterations of pre to apply
% last_time - an optional argument specifying the time in a time-varying
%    system to start applying pre.
%
% pre_union - the polyhedron that results from the iterated pre
%
%

if ~(target <= base_sys.Xterm)
    disp('Target set for pre union not contained in terminal set');
end

Kw_maps = cell(base_sys.T,1);
uc_maps = cell(base_sys.T,1);

n = base_sys.n;
m = base_sys.m;
l = base_sys.l;

% clear LP solver variables
yalmip('clear')
constraints = [];

x0 = sdpvar(n,1);
lambda = sdpvar(base_sys.T,1);
x0_c = sdpvar(n,base_sys.T);
constraints = [sum(lambda) == 1, lambda >= 0, x0 == sum(x0_c,2)];

for start_time = 1:base_sys.T
    
    horizon = base_sys.T - start_time + 1;
    sys = LTVSSys.constructTruncatedSystem(base_sys,start_time,horizon, target);
    
    initialConditions = computeInitialConditions(sys, target, sys.T, 0, 1);
    admissibleTrajectories = computeLiftedAdmissibleTrajectories(sys,sys.T);

    % We will search for control policies that are affine in the initial state
    % and the disturbance history given the switching history.
    % These containers map partial sequences to the appropriate decision
    % variable corresponding to the affine control law
    %       u = Kw * w + uc;
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
            uc = sdpvar(m,1,'full');
            uc_map(seq) = uc;
        end
    end
    Kw_maps{start_time} = Kw_map;
    uc_maps{start_time} = uc_map;


    % For each switching sequence of full length we add appropriate constraints
    total_sequences = sys.sequences{1,horizon+1};
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
            T*Hbar == Gu * Kw_inst + Gw * lambda(start_time);
            T*hbar <= lambda(start_time) * g - Gu * uc_inst - Gx * x0_c(:,start_time);
            T >= 0;
        ];
    end
end
% finally we solve the linear feasibility program with no objective
%diagnostics = optimize(constraints, [], yalmipOptions);

%(Polyhedron('H',[Hbar, hbar]) < Polyhedron('H',[value([Gx + Gu * Kx_inst, Gu * Kw_inst + Gw]), value(g - gu_offset)]))
%plot(Polyhedron(value([Gx + Gu * Kx_inst, Gu * Kw_inst + Gw]), value(g - gu_offset)))

%{
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
%}


end

