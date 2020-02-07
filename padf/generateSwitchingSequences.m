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
