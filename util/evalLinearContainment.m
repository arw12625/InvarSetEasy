function [ratios,containment] = evalLinearContainment(innerConstraint, innerVar, outerConstraint, outerVar, numberTests, origin)
%evalLinearContainment Evaluate the approximation or closeness of the sets
%specified by the inner and outer linear constraints. This is done by
%randomly generating rays rooted at the origin and evaluating the distance
%of their intersection with both of the constraint sets.
%   innerConstraint - YALMIP linear constraints for the point specified by the YALMIP variable innerVar
%   outerConstraint - YALMIP linear constraints for the point specified by the YALMIP variable outerVar
%   numberTests - the number of rays to evaluate
%   origin - the point at which the rays originate

dim = size(innerVar,1);

rays = randomUnitVec(dim, numberTests);

innerObjective = -rays * innerVar;
outerObjective = -rays * outerVar;

options = sdpsettings('solver', 'gurobi', 'verbose', 1);

innerDiag = optimize(innerConstraint, innerObjective, options);
innerDists = zeros(numberTests,1);
for i = 1:numberTests
    innerDists(i) = value(innerObjective(i), i);
end

outerDiag = optimize(outerConstraint, outerObjective, options);
outerDists = zeros(numberTests,1);
for i = 1:numberTests
    outerDists(i) = value(outerObjective(i), i);
end

ratios = zeros(numberTests,1);
for i = 1:numberTests
    offset = rays(i,:) * origin;
    ratios(i) = (-innerDists(i) + offset) / (-outerDists(i) + offset);
end

containment = 1;
if any(ratios > 1) || any(ratios < 0)
    containment = 0;
end

end

