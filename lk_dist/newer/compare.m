%{
% lane keeping
N = 10;
scale = 4.95;
Omega = Omega_lk
%}
N = 5;
scale = 1;
%{
consensus 3x3 grid
N=20;
scale = 8;
%}
N=20;
scale = 11;
%% linear program method

start_time = cputime;

[diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(scale * Omega, pslsys, N, yalmipOptions);

1 - diagnostics.problem

comp_time_lin = cputime - start_time

%% projection method

start_time = cputime;

s = scale * Omega;
for i = 1:N
    [i, size(s.A,1), cputime - start_time]
    s = polyLinPre(plsys, s);
end

(scale * Omega) <= s

comp_time_proj = cputime - start_time