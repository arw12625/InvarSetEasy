function [pre_s] = polyLinPre(plsys,s)
%polyLinPre Compute the pre in the linear system
    % plsys - polytopic linear system
    % s - set to compute pre of

Uv = plsys.U.V;

d = computeDisturbanceOffsets(plsys,s);

inputVertexPolys = Polyhedron(size(Uv,1));

for i = 1:size(Uv,1)
    distPoly = Polyhedron('A', s.A, 'b', s.b + d);
    inputVertexPolys(i) = invAffineMap(distPoly,plsys.A,plsys.B*Uv(i,:)');
end

pre_s = PolyUnion(inputVertexPolys).convexHull();
pre_s = intersect(pre_s, plsys.X);

end

