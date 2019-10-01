function [Gx, Gu, Gw, g] = computeSwitchedTargetConstraint(pslsys, A_seq, f_seq, target, switch_seq, N)
%   computeSwitchedTargetConstraint Generates a linear constraint on the
%   joint initial state, input sequence, and disturbance sequence space
%   that corresponds to trajectories that reach the target under the system
%   dynamics with the specified switching sequence.
%
%   pslsys - the polytopic switched linear system
%   A_seq - a container mapping that contains the products of A matrices
%               indexed by the corresponding string of switching signals
%   f_seq - a container mapping that contains the sums of f vectors
%               appropriately scaled by A matrices
%               indexed by the corresponding string of switching signals
%   target - the target polytope
%   switch_seq - the switching sequence for the dynamics
%   N - the total number of steps to be included in the constraint
%           Ie the constraint matrix is padded with zeros to correspond to
%           a total number of N steps
%
%   Gx - part of constraint matrix corresponding to x
%   Gu - part of constraint matrix corresponding to u
%   Gw - part of constraint matrix corresponding to w
%   g - constraint vector
%
%   The results are to be interpreted as a polytope or linear constraint by
%       P = {(x,u,w) | Gx * x + Gu * u + Gw * w <= g}
%

init_con_mat = A_seq(switch_seq);
split_seq = split(switch_seq, ',');

nsteps = size(split_seq, 1);

input_con_mat = cell(1,N);
dist_con_mat = cell(1,N);

partial_seq = '';
for i = nsteps:-1:1
    index = str2num(split_seq{i});
    input_con_mat{i} =  A_seq(partial_seq) * pslsys.B{index};
    dist_con_mat{i} = A_seq(partial_seq) * pslsys.E{index};
    if i == nsteps
        partial_seq = split_seq{i};
    else
        partial_seq = strcat(split_seq{i}, ',', partial_seq);
    end
end
for i = (nsteps+1):N
    input_con_mat{i} = zeros(size(pslsys.B{1}));
    dist_con_mat{i} = zeros(size(pslsys.E{1}));
end

input_con_mat = cell2mat(input_con_mat);
dist_con_mat = cell2mat(dist_con_mat);

Gx = target.A * init_con_mat;
Gu = target.A * input_con_mat;
Gw = target.A * dist_con_mat;
g = target.b - target.A * f_seq(switch_seq);

end