function [sequences, A_seq, f_seq] = computeSwitchedSequence(A,f,N)
%computeSwitchedSequence For each switching sequence of length at most N
%   compute the corresponding product of A matrices
%
%   A - cell array of system A matrices indexed by the switching modes 1:ns
%   f - cell array of system f vectors indexed by the switching modes 1:ns
%   N - the maximal length of switching signals considered
%
%   A_seq - a container mapping that contains the products of A matrices
%               indexed by the corresponding string of switching signals
%   f_seq - a container mapping that contains the sums of f vectors
%               appropriately scaled by A matrices
%               indexed by the corresponding string of switching signals
%
%   For example for the switching sequence s='1012' one would have that
%       A_seq(s) = A{2} * A{1} * A{0} * A{1}
%       f_seq(s) = f{2} + A{2} * f{1} + A{2} * A{1} * f{0} + A{2} * A{1} * A{0} * f{1} 
%
sequences = cell(N+1,1);
sequences{1} = {''};
n = size(A{1},1);
ns = size(A,1);
A_seq = containers.Map({''},{eye(n)});
f_seq = containers.Map({''},{zeros(n,1)});

for i = 1:N
    last_seq = sequences{i};
    new_seq = cell(ns * size(last_seq,1),1);
    delim='';
    if i ~= 1
        delim = ',';
    end
    for j = 1:size(last_seq,1)
        s = last_seq{j};
        for t = 1:ns
            new_str =  strcat(s,delim,num2str(t));
            new_seq{t + (j-1)*ns} = new_str;
            A_seq(new_str) = A{t} * A_seq(s);
            f_seq(new_str) = f{t} + A{t} * f_seq(s);
        end
    end
    sequences{i+1} = new_seq;
end

end
