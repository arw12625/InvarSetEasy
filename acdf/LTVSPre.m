function [s] = LTVSPre(sys, target, t)
% LTVSpre computes the pre (predecesor) of the target set at time t

s = Polyhedron(zeros(1,sys.n), 1);
for mode = 1:sys.ns
    W = sys.WSigmamap{t,mode};
    if isempty(W)
        continue;
    end
    s = s & linearPre(target,  sys.XUmap{t}, W, ...
                        sys.Amap{mode},  sys.Bmap{mode}, sys.Emap{mode}, ...
                        sys.fmap{mode});
end

end

function [pre_s] = linearPre(target, XU, W, A, B, E, f)
%linearPre Compute the pre in the linear system

d = computeDisturbanceOffsets(target, W, E);

As = [target.A; target.Ae; -target.Ae];
bs = [target.b; target.be; -target.be];

n = size(A,1);

H = [As * A, As * B, bs + d - As * f;
     XU.H];
He = XU.He;

lift_pre_s = Polyhedron('H', H, 'He', He);
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

TA = [target.A; target.Ae; -target.Ae];
%Tb = [target.be; target.be; -target.be];

%size(TA)
%size(E)
if sum(sum(abs(-TA * E))) == 0
    d = zeros(size(TA,1), 1);
    return;
end

w = sdpvar(size(E,2),1,'full');

objective = -TA * E * w;

constraints = [W.A * w <= W.b, W.Ae * w == W.be];

options = sdpsettings('solver', 'gurobi', 'verbose', 0);

diagnostics = optimize(constraints, objective, options);

d = zeros(size(TA,1),1);

for i = 1:size(TA,1)
    selectsolution(i);
    d(i) = value(objective(i));
end

end


