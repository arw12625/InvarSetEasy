function at = computeLiftedAdmissibleTrajectories(sys, horizon)
    % computeLiftedAdmissibleTrajectories defines the polyhedron
    % representing all lifted trajectories satisfying system constraints
    % A lifted trajectory consists of an initial state, an input sequence,
    % and a disturbance sequence over the horizon specified.
    %
    %   sys - an LTVSSys describing the system
    %   horizon - the horizon for the trajectories (must be at most the horizon of the system)
    %
    %   admissibleTrajectories is a structure representing the admissible
    %      trajectories of the system over the horizon. It consists of a
    %      map from switching sequences with length given by the horizon to
    %      polyhedra of admissible lifted trajectories. The switching
    %      sequence sigma corresponds to the polyhedron given by
    %           admis = { (x0,u,w) | at.Gxmap(sigma) * x0 + at.Gumap(sigma) * u + at.Gwmap(sigma) * w <= at.gmap(sigma) }

    
assert(horizon <= sys.T);
assert(horizon > 0);

n = sys.n;
m = sys.m;
l = sys.l;

% Each partial switching sequence corresponds to a target constraint either
% enforcing safety or reaching the target of Omega
% These constraints are of the form
%       Ax * x + Au * u + Aw * w <= b
% Equality constraints are transformed into inequality constraints
%
% The following variables define maps from partial sequences to the
% constraint matrices/vectors. They are initialized with the constraint
% corresponding to the empty string ''.
Ax_map = containers.Map('KeyType','char','ValueType','any');
Au_map = containers.Map('KeyType','char','ValueType','any');
Aw_map = containers.Map('KeyType','char','ValueType','any');
b_map = containers.Map('KeyType','char','ValueType','any');

% We define these constraints for all lengths of sequences from 1 to N
for t = 1:(horizon+1)
    len_sequences = sys.sequences{1,t};
    for j = 1:length(len_sequences)
        seq = len_sequences{j};
        if t == horizon+1
            Ax = [sys.Xterm.A; sys.Xterm.Ae; -sys.Xterm.Ae];
            Au = zeros(size(Ax,1),horizon*m);
            b = [sys.Xterm.b; sys.Xterm.be; -sys.Xterm.be];
        else
            Ax = [sys.XUmap{t}.A(:,1:n); sys.XUmap{t}.Ae(:,1:n); -sys.XUmap{t}.Ae(:,1:n)];
            Au = [zeros(size(Ax,1), (t-1)*m), ...
                 [sys.XUmap{t}.A(:,n+(1:m)); sys.XUmap{t}.Ae(:,n+(1:m)); -sys.XUmap{t}.Ae(:,n+(1:m))], ...
                  zeros(size(Ax,1), (horizon-t)*m)];
            b = [sys.XUmap{t}.b; sys.XUmap{t}.be; -sys.XUmap{t}.be];
        end
        
        Ax_dyn = Ax * sys.getSequenceA(seq);
        Au_dyn = cell(1,horizon);
        Aw_dyn = cell(1,horizon);
        
        modes = LTVSSys.getModesFromSequence(seq);
        for i = 1:length(modes)
            mode = modes(i);
            suffix = LTVSSys.getSequenceFromModes(modes(i+1:end));
            Au_dyn{i} = Ax * sys.getSequenceA(suffix) * sys.Bmap{mode};
            Aw_dyn{i} = Ax * sys.getSequenceA(suffix) * sys.Emap{mode};
        end
        for i = (length(modes)+1):horizon
            Au_dyn{i} = Ax * zeros(n,m);
            Aw_dyn{i} = Ax * zeros(n,l);
        end
        Au_dyn = Au + cell2mat(Au_dyn);
        Aw_dyn = cell2mat(Aw_dyn);
        
        Ax_map(seq) = Ax_dyn;
        Au_map(seq) = Au_dyn;
        Aw_map(seq) = Aw_dyn;
        b_map(seq) = b - Ax * sys.getSequencef(seq);
    end
end

% For each switching sequence of full length we add appropriate constraints
total_sequences = sys.sequences{1,horizon+1};
Ax_mapTotal = containers.Map('KeyType','char','ValueType','any');
Au_mapTotal = containers.Map('KeyType','char','ValueType','any');
Aw_mapTotal = containers.Map('KeyType','char','ValueType','any');
b_mapTotal = containers.Map('KeyType','char','ValueType','any');
for j = 1:length(total_sequences)

    sequence = total_sequences{j};
    
    % We build the constraint matrices in block form
    AxTotal = cell(horizon+1,1);
    AuTotal = cell(horizon+1,1);
    AwTotal = cell(horizon+1,1);
    bTotal = cell(horizon+1,1);
    
    modes = LTVSSys.getModesFromSequence(sequence);
    
    for t = 1:(horizon+1)
        prefix = LTVSSys.getSequenceFromModes(modes(1:(t-1)));
        AxTotal{t} = Ax_map(prefix);
        AuTotal{t} = Au_map(prefix);
        AwTotal{t} = Aw_map(prefix);
        bTotal{t} = b_map(prefix);
    end
    Ax_mapTotal(sequence) = cell2mat(AxTotal);
    Au_mapTotal(sequence) = cell2mat(AuTotal);
    Aw_mapTotal(sequence) = cell2mat(AwTotal);
    b_mapTotal(sequence) = cell2mat(bTotal);
end

at.Ax_map = Ax_mapTotal;
at.Au_map = Au_mapTotal;
at.Aw_map = Aw_mapTotal;
at.b_map = b_mapTotal;
at.horizon = horizon;

end