function [diagnostics] = computeController(sys,cfamily,uncertainty,admissible, yalmipOptions)
%COMPUTEAFFINECONTROLLER Summary of this function goes here
%   Detailed explanation goes here

horizon = cfamily.horizon;

A_un = [uncertainty.A; uncertainty.Ae; -uncertainty.Ae];
b_un = [uncertainty.b; uncertainty.be; -uncertainty.be];

A_ad = [admissible.A; admissible.Ae; -admissible.Ae];
b_ad = [admissible.b; admissible.be; -admissible.be];

sysmap = cfamily.computeSystemMap();

T = sdpvar(size(A_ad,1), size(A_un,1),'full');

F = [sysmap.Sx,sysmap.Sw,sysmap.Sv;
     sysmap.Ux,sysmap.Uw,sysmap.Uv];
f = [sysmap.s;
     sysmap.u];

constraints = [cfamily.constraints;
               T >= 0;
               T * A_un == A_ad * F;
               T * b_un <= b_ad - A_ad * f];

diagnostics = optimize(constraints, [], yalmipOptions);

end

