function iter_set = computeOuterInvar(sys,target,steps,last_time)
% computeIteratedPre computes the set that is the result of iteratively
% applying the pre operator to the target set.
% 
% If the target set is the safe set, this yields an outer approximation of
% the maximal invariant set.
%
% This methods inteprets the dynamics of the system as periodic.
%
% sys - the LTVSSys representing system dynamics
% target - the polyhedron to apply pre to
% steps - the number of iterations of pre to apply
% last_time - an optional argument specifying the time in a time-varying
%    system to start applying pre.
%
% iter_set - the polyhedron that results from the iterated pre
%
s = target;

if nargin == 3
    last_time = mod(steps, sys.T) + 1;
end

t = last_time;
for i = 1:steps
    s = s & sys.pre(s, t);
    t = t - 1;
    if t < 1
        t = sys.T;
    end
end

iter_set = s;

end

