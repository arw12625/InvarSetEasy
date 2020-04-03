function [isReach] = isAffineBackwardsReachable(sys,max_horizon,x0,XU,XF,W,V,options)
%COMPUTEAFFINECONTROLLER Summary of this function goes here
%   Detailed explanation goes here

if isEmptySet(W) && isEmptySet(V)
    error('computeControllerForPointNoDist cannot be use if the uncertainty set is empty');
end

x0_t = sdpvar(sys.n,max_horizon,'full');
lambda = sdpvar(max_horizon, 1);
constraints = [sum(lambda) == 1, lambda >= 0, x0 == sum(x0_t,2)];
cfamilies = cell(max_horizon,1);

for horizon = 1:max_horizon

    cfamily = ControlFamily(sys,horizon,true,true,false);
    cfamilies{horizon} = cfamily;

    uncertainty = polyhedronPower(W,horizon)*polyhedronPower(V,horizon);
    A_un = [uncertainty.A; uncertainty.Ae; -uncertainty.Ae];
    b_un = [uncertainty.b; uncertainty.be; -uncertainty.be];

    A_XU = [XU.A; XU.Ae; -XU.Ae];
    b_XU = [XU.b; XU.be; -XU.be];
    A_XU_x = A_XU(:,1:sys.n);
    A_XU_u = A_XU(:,sys.n+(1:sys.m));
    A_XF = [XF.A; XF.Ae; -XF.Ae];
    b_XF = [XF.b; XF.be; -XF.be];
    
    A_ad = [blkdiag(kron(eye(horizon),A_XU_x), A_XF), ...
            [kron(eye(horizon),A_XU_u); zeros(size(A_XF,1),horizon*sys.m)]];
    b_ad = [repmat(b_XU,horizon,1);b_XF];
    
    % w_tmat is a hack, sorry
    [sysmap, x_tmat, w_tmat] = cfamily.computeSystemMap();

    T = sdpvar(size(A_ad,1), size(A_un,1),'full');

    F = [sysmap.Sw + (-1+lambda(horizon))*w_tmat,sysmap.Sv;
         sysmap.Uw,sysmap.Uv];
    f = [sysmap.s + (sysmap.Sx + (-1+lambda(horizon))*x_tmat) * x0;
         sysmap.u];

    constraints = [constraints;
                   T >= 0;
                   T * A_un == A_ad * F;
                   T * b_un <= b_ad * lambda(horizon) - A_ad * f];

end

diagnostics = optimize(constraints,[],options);
isReach = diagnostics.problem == 0;

end

