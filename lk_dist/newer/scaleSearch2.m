
yalmip('clear')
clc

%% Load system variables


%% Compute seed set

%{
Napprox = 6;
%Omega = computeOuterApproxInvariantDist(X, U, W, Napprox, A, B, C);
%Omega = Polyhedron.unitBox(4);
load('working_omega', 'Omega');
% Need to compute the V-rep here, even though it is not used
% This removes redundancy in the H rep and simplifies it
% There is probably an MPT call meant for this instead
% I could replace the saved polytope with the fixed one, but this serves as
%   a reminder that the simplicity of the representation is crucial, and
%   that computing the V-rep recalculates the H-rep.
Omega.V;
%}

%% Search scalings
%{
4d lane keeping
N = 5;
scaleMult = 2.5;
scaleBase = 1.08;
scaleLogRange = 0:10;
%}

N = 5;
scaleMult = 1;
scaleBase = 1;
scaleLogRange = 0;

yalmipOptions = sdpsettings('solver', 'gurobi', 'verbose', 1);
success = [];
times = [];
for i = scaleLogRange
    starttime = cputime;
    scale = scaleMult * scaleBase^i;
    
    [diagnostics, sequences, Kx_map, Kw_map, uc_map] = evalSwitchedInvariance(scale * Omega_lk, pslsys, N, yalmipOptions);
    if diagnostics.problem == 0
        success = [success, scale];
    end
    times = [times, cputime - starttime];
end
success
times
