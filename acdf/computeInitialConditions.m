function initialConditions = computeInitialConditions(sys, initialStates, horizon, includeInitState, includeDisturbance)
    % computeInitialConditions determines initial conditions for a system
    % given by a map from switching sequences with length given by the horizon
    % to a polyhedron of initial states and possible disturbance sequences
    %
    %   sys - an LTVSSys representing the system
    %   initialStates - a polyhedron representing initial states
    %   horizon - the horizon for the initial conditions
    %   includeInitState - whether or not to include the initial state
    %   includeDisturbance - whether or not to include the disturbance
    %
    %   initialConditions is a structure mapping a switching sequence
    %   sigma to initial conditions by 
    %       initcond = { (x0,w) | H(sigma) * [x0;w] <= h(sigma) }
    %
    %
assert(horizon <= sys.T);
assert(horizon > 0);

if nargin <= 3
    includeInitState = 1;
    includeDisturbance = 1;
end
if nargin == 4
    includeDisturbance = 1;
end

initPolyMap = containers.Map('KeyType','char','ValueType','any');

total_sequences = sys.sequences{1,horizon+1};
for j = 1:size(total_sequences, 1)

    sequence = total_sequences{j};
    modes = LTVSSys.getModesFromSequence(sequence);
    
    if includeInitState
        A = initialStates.A;
        b = initialStates.b;
        Ae = initialStates.Ae;
        be = initialStates.be;
    else
        A = [];
        b = [];
        Ae = [];
        be = [];
    end
    if includeDisturbance
        for t=1:horizon
            mode = modes(t);
            if ~isEmptySet(sys.WSigmamap{t,mode})
                A = blkdiag(A,sys.WSigmamap{t,mode}.A);
                b = [b;sys.WSigmamap{t,mode}.b];
                Ae = blkdiag(Ae,sys.WSigmamap{t,mode}.Ae);
                be = [be;sys.WSigmamap{t,mode}.be];
            end
        end
    end
    
    initPolyMap(sequence) = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be);
end

initialConditions.horizon = horizon;
initialConditions.initPolyMap = initPolyMap;

end

