function [outerInvar] = computeOuterApproxInvariantDist(X,U,W,N,A,B,C)
t = X;
iter = N;

for i = 1:iter
    t = t & polyPreDist(A,B,C,X,U,W,t);
end

outerInvar = t;

end

