function [v] = randomUnitVec(dim, numVec)
%randomUnitVec Generates uniformly distributed unit vectors
%   dim - Dimension of vectors requested
%   numVec - Number of vectors requested

v = randn(numVec,dim);
v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));

end

