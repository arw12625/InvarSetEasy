function [pinv] = computeConvInvar(pslsys,Omega,N)

polys(1) = Polyhedron;
s = Omega;

for i = 1:N
    s = polySwitchedLinPre(pslsys, s);
    polys(i) = s;
end

pinv = convexHull(PolyUnion(polys));

end

