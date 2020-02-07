function [reachMap, totalReach] = computeControllerForwardReachableSet(sys, initSet, affineController)
    % reachableSetAffineController determines which states are reachable by
    % the given affineController from the initial state set.
    % This function returns a map (reachMap) from switching sequences to the states
    % reachable under this sequence.
    % The function also returns the convex hull of all reachable states (totalReach).
    
    horizon = affineController.horizon;
    
    delim = LTVSSys.getDelimiter();
    
    % map from switching sequence to polyhedron of possible current state,
    % initial state, and disturbance sequence tuples.
    aug_reach = containers.Map('KeyType','char','ValueType','any');
    
    % map from switching sequence to polyhedron of possible current states
    reachMap = containers.Map('KeyType','char','ValueType','any');
    
    % initialize maps with initial set
    aug_reach('') = [eye(sys.n);eye(sys.n)] * initSet;
    reachMap('') = initSet;
    
    totalReach = initSet;
    for t=1:horizon
        len_sequences = sys.sequences{1,t};
        for j = 1:length(len_sequences) % iterate over all sequences with length t
            seq = len_sequences{j};
            for mode = 1:sys.ns
                W = sys.WSigmamap{t,mode};
                if isempty(W)
                    continue;
                end
                % (x_t, x_0, w_0,..., w_{t-1}) -> (x_{t+1},x_0,w_0,..., w_{t-1}, w_t)
                %
                % x_{t+1} = sys.Amap{mode} * x_t + sys.Bmap{mode} *
                % (affineController.Kx_map(seq) * x_0 + 
                % affineController.Kw_map(seq) * (w_0,...,w_{t-1}) + 
                % affineController.uc_map(seq)) + 
                % sys.Emap{mode} * (w_0,...,w_{t-1}) + sys.fmap{mode}
                
                % Construct polyhedron matrices for augmented reachable set
                % under this sequence
                sW = aug_reach(seq) * W;
                pA = cell(3,3);
                pb = cell(3,1);
                pA{1,1} = sys.Amap{mode};
                pA{1,2} = sys.Bmap{mode} * affineController.Kx_map(seq);
                if t == 1
                    trunc_Kw = zeros(sys.m, 0);
                else
                    trunc_Kw = affineController.Kw_map(seq);
                    trunc_Kw = trunc_Kw(:,1:((t-1)*sys.l));
                end
                pA{1,3} = [sys.Bmap{mode} * trunc_Kw, sys.Emap{mode}];
                pb{1} = sys.fmap{mode} + sys.Bmap{mode} * affineController.uc_map(seq);
                
                pA{2,1} = zeros(sys.n);
                pA{2,2} = eye(sys.n);
                pA{2,3} = zeros(sys.n, t * sys.l);
                pb{2} = zeros(sys.n, 1);
                
                pA{3,1} = zeros(t * sys.l, sys.n);
                pA{3,2} = zeros(t * sys.l, sys.n);
                pA{3,3} = eye(t * sys.l);
                pb{3} = zeros(t * sys.l, 1);
                
                pA = cell2mat(pA);
                pb = cell2mat(pb);
                newS = (pA * sW) + pb;
                
                new_seq = strcat(seq, num2str(mode), delim);
                aug_reach(new_seq) = newS;
                reachMap(new_seq) = newS.projection(1:sys.n);
                totalReach(end+1) = reachMap(new_seq);
            end
        end
    end
    totalReach = PolyUnion(totalReach);
    totalReach = totalReach.convexHull();
end

