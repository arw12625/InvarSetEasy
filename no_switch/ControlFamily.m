classdef ControlFamily
    %CONTROLFAMILY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sys
        horizon
        parameters
        constraints
        allow_offset
        allow_feedback
        allow_feedback_init
    end
    
    methods
        function obj = ControlFamily(sys, horizon,allow_offset, allow_feedback, allow_feedback_init)
            %CONTROLFAMILY Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 2
                allow_offset = 1;
                allow_feedback = 1;
                allow_feedback_init = 1;
            end
            if nargin == 3
                allow_feedback = 1;
                allow_feedback_init = 1;
            end
            if nargin == 4
                allow_feedback_init = allow_feedback;
            end
            obj.sys = sys;
            obj.horizon = horizon;
            obj.allow_offset = allow_offset;
            obj.allow_feedback= allow_feedback;
            obj.allow_feedback_init = allow_feedback_init;
            if allow_feedback
                Kyt = sdpvar(sys.m * horizon, sys.k * horizon, 'full');
                for t = 1:horizon
                    Kyt(t+(0:sys.m-1), t*sys.k+1:sys.k*horizon) = zeros(sys.m, sys.k*(horizon - t));
                end
                if ~allow_feedback_init
                    Kyt(:,1:sys.k) = zeros(sys.m*horizon, sys.k);
                end
            else
                Kyt = zeros(sys.m * horizon, sys.k * horizon);
            end
            if allow_offset
                uc = sdpvar(sys.m * horizon, 1, 'full');
            else
                uc = zeros(sys.m * horizon, 1);
            end
            obj.parameters.Kyt = Kyt;
            obj.parameters.uc = uc;
            
            obj.constraints = [];
        end
        
        function [sysmap, x_tmat, w_tmat] = computeSystemMap(obj)
            [sysmap, x_tmat, w_tmat] = obj.sys.computeSystemMapModified(obj.horizon, obj.parameters.Kyt, obj.parameters.uc);
        end
        
        %{
        function [Sx,Sw,Sv,s,Ux,Uw,Uv,u] = computeStateInputAffineMap(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            sys1 = obj.sys;
            obs_tmat = LTISys.computeToeplitzMat(sys1.A,eye(sys1.n),sys1.C,obj.horizon);
            Uxw = obj.parameters.Kyt * [obs_tmat,zeros(obj.horizon * sys1.k,sys1.n)];
            Ux = Uxw(:,1:sys1.n);
            Uw = Uxw(:,sys1.n+1:end) * kron(eye(obj.horizon),sys1.E);
            Uv = obj.parameters.Kyt * kron(eye(obj.horizon),sys1.G);
            u = obj.parameters.uc;
            
            xw_tmat = LTISys.computeToeplitzMat(sys1.A,eye(sys1.n),eye(sys1.n),obj.horizon+1);
            x_tmat = xw_tmat(:,1:sys1.n);
            w_tmat = xw_tmat(:,sys1.n+1:end) * kron(eye(obj.horizon),sys1.E);
            u_tmat = [zeros(sys1.n,obj.horizon * sys1.m); LTISys.computeToeplitzMat(sys1.A,sys1.B,eye(sys1.n),obj.horizon)];
            Sx = x_tmat + u_tmat * Ux;
            Sw = w_tmat + u_tmat * Uw;
            Sv = u_tmat * Uv;
            s = u_tmat * u;
        end
        %}
        function [Ky,uc] = computeOutputFeedbacks(obj)
            [Ky,uc] = obj.sys.computeOutputFeedbacks(obj.horizon, value(obj.parameters.Kyt),value(obj.parameters.uc));
        end
        
    end
    
end

