function [outerInvar] = computeOuterApproxInvariant(X,U,N,A,B)

%{
t = X;
iter = N;

for i = 1:iter
    t = t & polyPre(A,B,X,U,t);
end

outerInvar = t;
%}

system = LTISystem('A',A,'B',B);

outerInvar = system.invariantSet('X',X,'U',U,'maxIterations',N);

end

