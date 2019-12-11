function initialConditions = computeInitialConditions(sys, initialStates, horizon)
    % computeInitialConditions determines initial conditions for a system
    % given by a map from switching sequences with length given by the horizon
    % to a polyhedron of initial states and possible disturbance sequences
    %
    %   sys - an LTVSSys representing the system
    %   initialStates - a polyhedron representing initial states
    %   horizon - the horizon for the initial conditions
    %
    %   initialConditions is a structure mapping a switching sequence
    %   sigma to initial conditions by 
    %       initcond = { (x0,w) | H(sigma) * [x0;w] <= h(sigma) }
    %
    %
assert(horizon <= sys.T);
assert(horizon > 0);

Hmap = containers.Map('KeyType','char','ValueType','any');
hmap = containers.Map('KeyType','char','ValueType','any');

total_sequences = sys.sequences{1,horizon+1};
for j = 1:size(total_sequences, 1)

    sequence = total_sequences{j};
    modes = LTVSSys.getModesFromSequence(sequence);
    
    H = initialStates.A;
    h = initialStates.b;
    for t=1:horizon
        mode = modes(t);
        if ~isEmptySet(sys.WSigmamap{t,mode})
            H = blkdiag(H,sys.WSigmamap{t,mode}.A);
            h = [h;sys.WSigmamap{t,mode}.b];
        end
    end
    
    Hmap(sequence) = H;
    hmap(sequence) = h;
end

initialConditions.Hmap = Hmap;
initialConditions.hmap = hmap;
initialConditions.horizon = horizon;
end

