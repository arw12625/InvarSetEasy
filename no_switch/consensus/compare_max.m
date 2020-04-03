ssys = LTVSSys.constructLTISys(T,A,B,E,f,X*U,X,W);
tic
invar = computeOuterInvar(ssys,X,T);
toc

%%
XF = Polyhedron('A',P,'b',q);
isReachArray = zeros(size(invar.V,1),1);
for i = 1:size(invar.V,1)
    x0_test = invar.V(i,:)';
    isReach = isAffineBackwardsReachable(sys,15,x0_test,X*U,XF,W,V,options)
    isReachArray(i) = isReach;
end