%% System definitions

% This script demonstrates various invariant set computations for a 
% discrete time double integrator using this library. 

% System matrices
% x - state
% u - input
% w - disturbance
% x(t+1) = A x(t) + B u(t) + E w(t) + f
A = [1,1;0,1];
B = [0;1];
E = [0;1];
f = [0;0];

% System constraints as polyhedra
Xsafe = (1 * Polyhedron.unitBox(1)) * (1 * Polyhedron.unitBox(1));
Usafe = Polyhedron.unitBox(1);
XU = Xsafe * Usafe;
Xterm = Xsafe;
W = 0.25 * Polyhedron.unitBox(1);

% The horizon considered
T = 6;

% Construct a linear time-varying switched system (LTVSS) representing the linear
% time-invariant double integrator system.
% For now all methods operate on LTVSS to avoid redundant functionality
sys = LTVSSys.constructLTISys(T,A,B,E,f,XU,Xterm,W);

%% Invariant set computations

% The seed set used for computation
Omega = (0.4* Polyhedron.unitBox(1)) * (0.4 * Polyhedron.unitBox(1));

% Determine if there is an affine disturbance feedback controller that can
% drive the system starting in Omega into Omega in exactly T steps subject
% to all disturbances and constraints.
% If the system is recurrent (satisfies this condition), then the method
% returns such a controller.
[is_rec, controller] = testAffineRecurrence(sys,Omega,T);

%%

% Determine all possible reachable states using this controller over the
% horizon. If the system is recurrent, then total reach set is invariant.
[reachMap, totalReach] = computeControllerForwardReachableSet(sys, Omega, controller);

%%

% Determine the states that can reach Omega within the horizon.
% (Note technically this set computes the convex hull of the union 
% of the sets pre^i(Omega) for i = 1:T which is slightly different)
% If the system is recurrent then this set is invariant.
preUnion = computeBackwardsReachableSet(sys, Omega, T);

%%

% Compute an outer approximation of the maximal invariant set by
% iteratively applying pre to the safe set.
%tic
outerInvar = computeOuterInvar(sys, Xsafe, 3);
%toc


%%

% Plot all of the sets computed
figure(T)
%plot([Xsafe, outerInvar, preUnion, totalReach, Omega])
%legend(["Safe Set", "Maximal Invariant", "Union of Pre", "Reach of Controller", "Seed Set"])
plot([Xsafe, outerInvar, preUnion, totalReach, Omega])
legend(["Safe Set X_s", "Maximal Invariant", "Union of Pre \Omega_N", "Reach of Controller", "Seed Set \Omega"])
title("Linear System Invariant Sets with N=6")
%{
% Plot the reachable sets
for t = 1:T
    figure(t)
    seqs = sys.sequences{1,t};
    plot(reachMap(seqs{1}))
end
%}

%%

[constraints, x0,lambda,x0_c,Kw_maps,uc_maps] = computeAffineBackwardsReachableSet(sys, Omega);
yalmipOptions = sdpsettings('verbose', 0, 'solver', ''); % options for the LP solver

infeas = 0;
for i = 1:size(outerInvar.V,2)
    mod_con = [constraints, x0 == outerInvar.V(i,:)'];
    diag = optimize(mod_con, [], yalmipOptions);
    if diag.problem ~= 0
        infeas = 1;
        break
    end
end

if infeas
    disp('Affine backwards reachable set is not maximal invariant set')
else
    disp('Affine backwards reachable set is maximal invariant set')
end