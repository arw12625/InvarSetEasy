function [pre_s] = polyLinPre(plsys,s)
%polyLinPre Compute the pre in the linear system
    % plsys - polytopic linear system
    % s - set to compute pre of

d = computeDisturbanceOffsets(plsys,s);

Ax = plsys.X.A;
bx = plsys.X.b;
As = s.A;
bs = s.b;
Au = plsys.U.A;
bu = plsys.U.b;

n = size(Ax,2);
m = size(Au,2);

H = [As * plsys.A, As * plsys.B, bs + d;
     Ax, zeros(size(Ax,1), m), bx; 
     zeros(size(Au,1), n), Au, bu];

lift_pre_s = Polyhedron('H', H);
pre_s = lift_pre_s.projection(1:n);

end

