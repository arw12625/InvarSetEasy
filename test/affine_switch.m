N = 1;
n = 1;
m = 1;
l = 1;
ns = 2;

A = cell(ns,1);
A{1} = 1;
A{2} = 1;

B = cell(ns,1);
B{1} = 1;
B{2} = 1;

E = cell(ns,2);
E{1} = 0;
E{2} = 0;

f = cell(ns,1);
f{1} = 2;
f{2} = 1;

Omega = 1 * Polyhedron.unitBox(n);
X = 10 * Polyhedron.unitBox(n);
U = 2 * Polyhedron.unitBox(m);
W = 0 * Polyhedron.unitBox(l);

pslsys = PolySwitchLinSys(A,X,B,U,E,W,f);
yalmipOptions = sdpsettings('verbose', 1);

[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(Omega, pslsys, N, yalmipOptions);

diagnostics.problem
%{
if diagnostics.problem == 0
    for i = 1:(size(sequences, 1) - 1)
        len_seq = sequences{i};
        for j = 1:size(len_seq, 1)
            seq = len_seq{j}
            value(Kx_map(seq))
            value(Kw_map(seq))
        end
    end
end
%}