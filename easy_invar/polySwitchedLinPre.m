function [pre_s] = polySwitchedLinPre(pslsys,s)
%polyLinPre Compute the pre in the switched linear system
    % pslsys - polytopic switched linear system
    % s - set to compute pre of

poly = s;

%compute pre under each mode and intersect
for i = 1:pslsys.ns
    mode_sys = PolyLinSys( ...
                    pslsys.A{i}, ...
                    pslsys.X, ...
                    pslsys.B{i}, ...
                    pslsys.U, ...
                    pslsys.E{i}, ...
                    pslsys.W, ...
                    pslsys.f{i} ...
               );
           
    pre = polyLinPre(mode_sys, s);
    poly = intersect(poly, pre);

end

pre_s = poly;

end

