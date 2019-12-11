function [s] = LTVSPre(sys, target, t)
% LTVSpre computes the pre (predecesor) of the target set at time t

s = Polyhedron(zeros(1,sys.n), 1);
for mode = 1:sys.ns
    W = sys.WSigmamap{t,mode};
    if isempty(W)
        continue;
    end
    s = s & regularpre(target,  sys.XUmap{t}, W, ...
                        sys.Amap{mode},  sys.Bmap{mode}, sys.Emap{mode}, ...
                        sys.fmap{mode});
end

end

function [pre_s] = regularpre(target, XU, W, A, B, E, f)
%regularpre Compute the pre in the linear system

d = computeDisturbanceOffsets(target, W, E);

As = target.A;
bs = target.b;

n = size(A,1);

H = [As * A, As * B, bs + d - As * f;
     XU.H];

lift_pre_s = Polyhedron('H', H);
pre_s = lift_pre_s.projection(1:n);

end


function [d] = computeDisturbanceOffsets(target, W, E)
% computeDisturbanceOffsets Computes the halfspace offsets that represent
% the effect of the disturbance upon Pre operation
%   target - polyhedron to compute pre of
%   W - polyhedron of disturbances
%   
%   dynamics given by x(t+1) = A * x(t) + B * u(t) + E * w(t) + f
%
%   d - disturbance offset

TA = target.A;
WA = W.A;
Wb = W.b;

size(TA)
size(E)
if sum(sum(abs(-TA * E))) == 0
    d = zeros(size(TA,1), 1);
    return;
end

w = sdpvar(size(E,2),1,'full');

objective = -TA * E * w;

constraints = (WA * w <= Wb);

options = sdpsettings('solver', 'gurobi', 'verbose', 0);

diagnostics = optimize(constraints, objective, options);

d = zeros(size(TA,1),1);

for i = 1:size(TA,1)
    selectsolution(i);
    d(i) = value(objective(i));
end

end


