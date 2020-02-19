classdef ControlFamily
    %CONTROLFAMILY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        horizon
        parameters
        constraints
    end
    
    methods
        function obj = ControlFamily(sys, horizon)
            %CONTROLFAMILY Construct an instance of this class
            %   Detailed explanation goes here
            obj.horizon = horizon;
            Kyt = sdpvar(sys.m * horizon, sys.k * horizon, 'full');
            uc = sdpvar(sys.m * horizon, 1, 'full');
            for t = 1:horizon
                Kyt(t+(0:sys.m-1), t*sys.k+1:sys.k*horizon) = zeros(sys.m, sys.k*(horizon - t));
            end
            obj.parameters.Kyt = Kyt;
            obj.parameters.uc = uc;
            
            %obj.constraints = [uc == 0 * uc];
            obj.constraints = [];
        end
        
        function [Sx,Sw,Sv,s,Ux,Uw,Uv,u] = computeStateInputAffineMap(obj,sys)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obs_tmat = LTISys.computeToeplitzMat(sys.A,eye(sys.n),sys.C,obj.horizon);
            Uxw = obj.parameters.Kyt * [obs_tmat,zeros(obj.horizon * sys.k,sys.n)];
            Ux = Uxw(:,1:sys.n);
            Uw = Uxw(:,sys.n+1:end);
            Uv = obj.parameters.Kyt;
            u = obj.parameters.uc;
            
            xw_tmat = LTISys.computeToeplitzMat(sys.A,eye(sys.n),eye(sys.n),obj.horizon+1);
            x_tmat = xw_tmat(:,1:sys.n);
            w_tmat = xw_tmat(:,sys.n+1:end);
            u_tmat = [zeros(sys.n,obj.horizon * sys.m); LTISys.computeToeplitzMat(sys.A,sys.B,eye(sys.n),obj.horizon)];
            Sx = x_tmat + u_tmat * Ux;
            Sw = w_tmat + u_tmat * Uw;
            Sv = u_tmat * Uv;
            s = u_tmat * u;
        end
        
    end
end

