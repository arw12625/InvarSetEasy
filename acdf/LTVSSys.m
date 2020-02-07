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
                generateSwitchingSequences(obj.T, obj.n, obj.ns, ...
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

