function [diagnostics] = computeControllerForPointNoDist(sys,cfamily,x0,admissible, yalmipOptions)
%COMPUTEAFFINECONTROLLER Summary of this function goes here
%   Detailed explanation goes here

horizon = cfamily.horizon;

A_ad = [admissible.A; admissible.Ae; -admissible.Ae];
b_ad = [admissible.b; admissible.be; -admissible.be];

sysmap = cfamily.computeSystemMap();

F = [sysmap.Sx
     sysmap.Ux];
f = [sysmap.s;
     sysmap.u];

constraints = [cfamily.constraints;
               A_ad * F * x0 <= b_ad - A_ad * f];

diagnostics = optimize(constraints, [], yalmipOptions);

end

