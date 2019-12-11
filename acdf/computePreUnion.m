function [pre_union] = computePreUnion(sys, target, steps, last_time)
% computePreUnion computes the convex hull of the union of the sets produced by iteratively
% applying the pre operator to the target set for the given number of steps
%
%
% sys - the LTVSSys representing system dynamics
% target - the polyhedron to apply pre to
% steps - the number of iterations of pre to apply
% last_time - an optional argument specifying the time in a time-varying
%    system to start applying pre.
%
% pre_union - the polyhedron that results from the iterated pre
%
%
s = target;
preSets = Polyhedron();

if nargin == 3
    last_time = sys.T;
end

t = last_time;
for i = 1:steps
    s = LTVSPre(sys, s, t);
    preSets(end+1) = s;
    t = t - 1;
    if t < 1
        t = sys.T;
    end
end

pre_union = PolyUnion(preSets);
pre_union = pre_union.convexHull();

end

