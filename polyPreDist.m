function [pre_s] = polyPreDist(A, B, C, X, U, W, s)
%polyPreDist computes the Pre or backwards step of a polytopic region in a linear system
%   A - transition matrix
%   B - input matrix
%   X - system domain polytope
%   W - disturbance polytope
%   U - input polytope 
%   s - polytope to compute Pre of

Uv = U.V;

d = computeDisturbanceOffsets(s,W,1,A,C);

inputVertexPolys = Polyhedron(size(Uv,1));

[s.b, d]

for i = 1:size(Uv,1)
    distPoly = Polyhedron('A', s.A, 'b', s.b + d);
    inputVertexPolys(i) = invAffineMap(distPoly,A,B*Uv(i,:)');
end

pre_s = PolyUnion(inputVertexPolys).convexHull();
pre_s = intersect(pre_s, X);

end

