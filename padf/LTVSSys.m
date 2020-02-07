classdef LTVSSys < matlab.mixin.Copyable
    % LTVSSys represents a linear time-varying switched system in discrete time
    %
    % More specifically, LTVSSys encompasses time-invariant affine switched
    % dynamics over a fixed horizon T with fixed time-varying joint
    % disturbance/switching mode constraints and joint state/input
    % constraints as well as a terminal state constraint.
    %
    % t < T
    % x(t+1) = Amap{sigma(t)} x(t) + Bmap{sigma(t)} u(t) + Emap{sigma(t)} w(t) + fmap{sigma(t)}
    % x(t) in R^n, u(t) in R^m, w(t) in R^l
    % sigma(t) in {1,2,...,ns}
    % (x(t),u(t)) in XUmap{t}
    % w(t) in WSigmamap{t,sigma(t)}
    % x(T) in Xterm
    %
    % Switching modes are represented by numbers 1:ns while switching
    % sequences are represented by strings of switching modes each followed
    % by a delimiter. For example the sequence corresponding to (1,2,4,1) 
    % would be "1,2,4,1,". The empty sequence () corresponds to the empty
    % string "".
    %
    % This class enumerates all possible switching sequences over the
    % horizon and the associated transition function matrices. These
    % quantities are indexed by the sequence. For instance for a switching
    % sequence of (1,2,4,1) the corresponding A transition matrix would be
    % Aseq("1,2,4,1,") = Amap{1} * Amap{4} * Amap{2} * Amap{1}
    %
    % Constraints are represented by MPT3 Polyhedrons
    %
    % An LTI system can be constructed by using the function constructLTISys
    
    properties (SetAccess = immutable)
        T % horizon
        Amap % cell array of A matrices indexed by switching mode
        Bmap % cell array of B matrices indexed by switching mode
        Emap % cell array of E matrices indexed by switching mode
        fmap % cell array of f vectors indexed by switching mode
        n % state dimension
        m % input dimension
        l % disturbance dimension
        ns % number of switching modes
        WSigmamap % cell array of disturbance constraint polyhedrons indexed by (time, switching mode)
        sequences % switching mode sequences
    end
    properties
        XUmap % cell array of joint state input constraint polyhedrons indexed by time
        Xterm % Terminal state constraint polyhedron
    end
    properties (SetAccess = private)
        Aseq % transition A matrices indexed by sequences
        fseq % transition f vectors indexed by sequences
    end
    methods
        function obj = LTVSSys(T,Amap,Bmap,Emap,fmap,XUmap,Xterm,WSigmamap)
            %PolySwitchLinSys Construct an instance of this class
            obj.T = T;
            obj.Amap = Amap;
            obj.Bmap = Bmap;
            obj.Emap = Emap;
            obj.fmap = fmap;
            
            obj.ns = length(Amap(1,:));
            
            % find possible first mode
            pos_mode = 0;
            for i = 1:obj.ns
                if ~isempty(Amap{1,i})
                    pos_mode = i;
                    break
                end
            end
            
            obj.n = size(Amap{pos_mode},2);
            obj.m = size(Bmap{pos_mode},2);
            obj.l = size(Emap{pos_mode},2);
            
            obj.XUmap = XUmap;
            obj.Xterm = Xterm;
            obj.WSigmamap = WSigmamap;
            
            obj.assertValidSystem();
            
            [obj.sequences, obj.Aseq, obj.fseq] = ...
                LTVSSys.generateSwitchingSequences(obj.T, obj.n, obj.ns, ...
                obj.Amap, obj.fmap, obj.WSigmamap);
        end
        
        function Am = getSequenceA(obj, seq)
            Am = obj.Aseq(seq);
        end
        
        function fv = getSequencef(obj, seq)
            fv = obj.fseq(seq);
        end
        
        function [] = assertValidSystem(obj)
            assert(obj.T > 0);
            assert(obj.ns > 0);
            assert(all(size(obj.Amap) == [obj.ns,1]));
            assert(all(size(obj.Bmap) == [obj.ns,1]));
            assert(all(size(obj.Emap) == [obj.ns,1]));
            assert(all(size(obj.fmap) == [obj.ns,1]));
            for mode = 1:obj.ns
                assert(all(size(obj.Amap{mode}) == [obj.n,obj.n]));
                assert(all(size(obj.Bmap{mode}) == [obj.n,obj.m]));
                assert(all(size(obj.Emap{mode}) == [obj.n,obj.l]));
                assert(all(size(obj.fmap{mode}) == [obj.n,1]));
            end
            assert(all(size(obj.XUmap) == [obj.T,1]));
            assert(all(size(obj.WSigmamap) == [obj.T,obj.ns]));
            for t = 1:obj.T
                assert(size(obj.XUmap{t}.A,2) == obj.n+obj.m);
                nodist = 1;
                for mode = 1:obj.ns
                    W = obj.WSigmamap{t,mode};
                    if ~isempty(W)
                        nodist = 0;
                        assert(size(W.A,2) == obj.l);
                    end
                end
                assert(nodist == 0);
            end
            assert(size(obj.Xterm.A,2) == obj.n);
        end
        
        function [s] = pre(obj, target, t)
            % pre computes the pre (predecesor) of the target set at time t

            s = Polyhedron(zeros(1,obj.n), 1);
            for mode = 1:obj.ns
                W = obj.WSigmamap{t,mode};
                if isempty(W)
                    continue;
                end
                s = s & LTVSSys.linearPre(target,  obj.XUmap{t}, W, ...
                                    obj.Amap{mode},  obj.Bmap{mode}, obj.Emap{mode}, ...
                                    obj.fmap{mode});
            end
        end
    end
    
    methods(Static, Access = private)

        function [sequences, A_seq, f_seq] = generateSwitchingSequences(T,n,ns, Amap, fmap, WSigmamap)
            % generateSwitchingSequences generates all possible switching sequences
            % and associated transition matrices for a system.
            %   T - horizon
            %   n - number of states
            %   ns - number of switching modes
            %   Amap - cell array of A matrices indexed by switching mode
            %   fmap - cell array of f matrices indexed by switching mode
            %   WSigmamap - cell array of disturbance polyhedrons indexed by time
            %      and switching mode
            %   
            %   sequences - cell array of list of switching sequence strings
            %      indexed by length.

            sequences = cell(T+1,T+1);
            A_seq = containers.Map({''},{eye(n)});
            f_seq = containers.Map({''},{zeros(n,1)});

            for st = 1:T
                sequences{st,1} = {''};
                for lt = 1:(T-st+1)
                    prev_sequences = sequences{st,lt};
                    new_sequences = cell(0,0);
                    delim = LTVSSys.getDelimiter();
                    for j = 1:size(prev_sequences,1)
                        seq = prev_sequences{j};
                        for mode = 1:ns
                            if ~isempty(WSigmamap{st+lt-1,mode}) && ~isEmptySet(WSigmamap{st+lt-1,mode})
                                new_str =  strcat(seq,num2str(mode),delim);
                                new_sequences = {new_sequences; new_str};
                                if ~isKey(A_seq, new_str)
                                    A_seq(new_str) = Amap{mode} * A_seq(seq);
                                    f_seq(new_str) = fmap{mode} + Amap{mode} * f_seq(seq);
                                end
                            end
                        end
                    end
                    sequences{st,lt+1} = new_sequences(2:end);
                end
            end

        end
        
        function [pre_s] = linearPre(target, XU, W, A, B, E, f)
            %linearPre Compute the pre in the linear system

            d = LTVSSys.computeDisturbanceOffsets(target, W, E);
            As = [target.A; target.Ae; -target.Ae];
            bs = [target.b; target.be; -target.be];
            n = size(A,1);
            H = [As * A, As * B, bs + d - As * f;
                 XU.H];
            He = XU.He;

            lift_pre_s = Polyhedron('H', H, 'He', He);
            pre_s = lift_pre_s.projection(1:n);

        end


        function [d] = computeDisturbanceOffsets(target, W, E)
            % computeDisturbanceOffsets Computes the halfspace offsets that represent
            % the effect of the disturbance upon Pre operation
            %   target - polyhedron to compute pre of
            %   W - polyhedron of disturbances
            %   
            %   dynamics given by x(t+1) = A * x(t) + B * u(t) + E * w(t) + f
            %
            %   d - disturbance offset

            TA = [target.A; target.Ae; -target.Ae];
            %Tb = [target.be; target.be; -target.be];

            %size(TA)
            %size(E)
            if sum(sum(abs(-TA * E))) == 0
                d = zeros(size(TA,1), 1);
                return;
            end

            w = sdpvar(size(E,2),1,'full');
            objective = -TA * E * w;
            constraints = [W.A * w <= W.b, W.Ae * w == W.be];
            options = sdpsettings('solver', 'gurobi', 'verbose', 0);
            diagnostics = optimize(constraints, objective, options);

            d = zeros(size(TA,1),1);
            for i = 1:size(TA,1)
                selectsolution(i);
                d(i) = value(objective(i));
            end
        end
    end
    
    methods(Static)
        function delim = getDelimiter()
            %return the deilimiter used to separate sequences
            delim = ',';
        end
        
        function mode = getLastMode(s)
            %return the last mode of a sequence
            mustBeNonempty(s)
            modes = split(s,LTVSsys.delimiter);
            mode = modes(end);
        end
        
        function sequence = getSequenceFromModes(modes)
            %transform an array of modes into sequence string
            if isempty(modes)
                sequence = "";
            else
                sequence = sprintf("%d,",modes);
            end
        end
        
        function modes = getModesFromSequence(seq)
            %transforma a sequence string into an array of modes
            modes = str2num(seq);
        end
        
        function len = getSequenceLength(seq)
            %get the number of modes in a sequence
            len = length(str2num(seq));
        end
        
        function sys = constructLTISys(T,A,B,E,f,XU,Xterm,W)
            % construct an LTVSSys representing an LTI system
            %   A,B,E,f - matrices
            %   XU - a single time-invariant state-input constraint polyhedron
            %           or cell array of such polyhedrons indexed by time
            %   W - disturbance polyhedron
            ns = 1;
            Amap = cell(ns,1);
            Amap(:) = {A};
            Bmap = cell(ns,1);
            Bmap(:) = {B};
            Emap = cell(ns,1);
            Emap(:) = {E};
            fmap = cell(ns,1);
            fmap(:) = {f};
            
            XUmap = cell(T,1);
            if iscell(XU)
                XUmap = XU;
            else
                XUmap(:) = {XU};
            end
            
            WSigmamap = cell(T,ns);
            WSigmamap(:) = {W};
            
            sys = LTVSSys(T,Amap,Bmap,Emap,fmap,XUmap,Xterm,WSigmamap);

        end
        
        function seq_sys = constructLTVSysFromSeq(sys, seq, start_time, Xterm)
            
            T = LTVSSys.getSequenceLength(seq);
            modes = LTVSSys.getModesFromSequence(seq);
            
            WSigmamap = cell(T,sys.ns);
            for i = 1:T
                WSigmamap{i,modes(i)} = sys.WSigmamap{start_time+i,modes(i)};
            end
            
            seq_sys = LTVSSys(T,sys.Amap,sys.Bmap,sys.Emap,sys.fmap,sys.XUmap(start_time:start_time+T-1),Xterm,WSigmamap);
            
        end
        
        function trunc_sys = constructTruncatedSystem(sys,start_time,T, Xterm)
            
            time_range = start_time+(1:T)-1;
            trunc_sys = LTVSSys(T,sys.Amap,sys.Bmap,sys.Emap,sys.fmap,sys.XUmap(time_range),Xterm,sys.WSigmamap(time_range));
            
        end
        
    end
end

