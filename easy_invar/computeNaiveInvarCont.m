function [pinv] = computeNaiveInvarCont(pslsys, Omega, N, sequences, Kx_map, Kw_map, uc_map)

polys(1) = Omega;

for seq = sequences{N+1}
    s = Omega; 
    %compute invar set induced by controller
end

pinv = convexHull(PolyUnion(polys));


end

