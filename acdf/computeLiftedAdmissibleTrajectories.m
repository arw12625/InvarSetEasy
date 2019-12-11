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
%       Gx * x + Gu * u + Gw * w <= g
%
% The following variables define maps from partial sequences to the
% constraint matrices/vectors. They are initialized with the constraint
% corresponding to the empty string ''.
target_G_map_x = containers.Map('KeyType','char','ValueType','any');
target_G_map_u = containers.Map('KeyType','char','ValueType','any');
target_G_map_w = containers.Map('KeyType','char','ValueType','any');
target_g_map = containers.Map('KeyType','char','ValueType','any');

% We define these constraints for all lengths of sequences from 1 to N
for t = 1:(horizon+1)
    len_sequences = sys.sequences{1,t};
    for j = length(len_sequences)
        seq = len_sequences{j};

        if t == horizon+1
            Gx = sys.Xterm.A;
            Gu = zeros(size(Gx,1),horizon*m);
            g = sys.Xterm.b;
        else
            Gx = sys.XUmap{t}.A(:,1:n);
            Gu = [zeros(size(Gx,1), (t-1)*m),  sys.XUmap{t}.A(:,n+(1:m)), zeros(size(Gx,1), (horizon-t)*m)];
            g = sys.XUmap{t}.b;
        end
        
        Gbx = Gx * sys.getSequenceA(seq);

        Gbu = cell(1,horizon);
        Gbw = cell(1,horizon);

        modes = LTVSSys.getModesFromSequence(seq);
        for i = 1:length(modes)
            mode = modes(i);
            suffix = LTVSSys.getSequenceFromModes(modes(i+1:end));
            Gbu{i} = Gx * sys.getSequenceA(suffix) * sys.Bmap{mode};
            Gbw{i} = Gx * sys.getSequenceA(suffix) * sys.Emap{mode};
        end
        for i = (length(modes)+1):horizon
            Gbu{i} = Gx * zeros(n,m);
            Gbw{i} = Gx * zeros(n,l);
        end
        Gbu = Gu + cell2mat(Gbu);
        Gbw = cell2mat(Gbw);
        
        target_G_map_x(seq) = Gbx;
        target_G_map_u(seq) = Gbu;
        target_G_map_w(seq) = Gbw;
        target_g_map(seq) = g - Gx * sys.getSequencef(seq);
    end
end

% For each switching sequence of full length we add appropriate constraints
total_sequences = sys.sequences{1,horizon+1};
Gxmap = containers.Map('KeyType','char','ValueType','any');
Gumap = containers.Map('KeyType','char','ValueType','any');
Gwmap = containers.Map('KeyType','char','ValueType','any');
gmap = containers.Map('KeyType','char','ValueType','any');
for j = 1:size(total_sequences, 1)

    sequence = total_sequences{j};
    
    % We build the constraint matrices in block form
    target_G_x = cell(horizon+1,1);
    target_G_u = cell(horizon+1,1);
    target_G_w = cell(horizon+1,1);
    target_g = cell(horizon+1,1);
    
    modes = LTVSSys.getModesFromSequence(sequence);
    
    for t = 1:(horizon+1)
        prefix = LTVSSys.getSequenceFromModes(modes(1:(t-1)));
        target_G_x{t} = target_G_map_x(prefix);
        target_G_u{t} = target_G_map_u(prefix);
        target_G_w{t} = target_G_map_w(prefix);
        target_g{t} = target_g_map(prefix);
    end
    Gxmap(sequence) = cell2mat(target_G_x);
    Gumap(sequence) = cell2mat(target_G_u);
    Gwmap(sequence) = cell2mat(target_G_w);
    gmap(sequence) = cell2mat(target_g);
    
end

at.Gxmap = Gxmap;
at.Gumap = Gumap;
at.Gwmap = Gwmap;
at.gmap = gmap;
at.horizon = horizon;

end