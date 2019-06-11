clear all
yalmip clear
clc

%% Load system variables

load('A', 'A');
load('B', 'B');
load('E', 'E');
load('U', 'U');
load('X', 'X');
load('W', 'W');
plsys = PolyLinSys(A,X,B,U,E,W);

%% Compute seed set

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


%% Search scalings

%scale 0.4786 works well for N=6
%scale 0.71 works well for N=30

N = 6;
scaleMult = 0.47;
scaleBase = 1.0;
scaleLogRange = 0:0;

yalmipOptions = sdpsettings('solver', 'gurobi', 'verbose', 1);
success = [];
for i = scaleLogRange
    scale = scaleMult * scaleBase^i;
    
    [isInvar, diagnostics, Kx,Kw,u0] = evalEasyControlInvariance(scale * Omega, plsys, N, 1,0, yalmipOptions);
    if isInvar
        success = [success, scale];
    end
end
success

