classdef LTISys < matlab.mixin.Copyable
    
    properties (SetAccess = immutable)
        A % A matrix
        B % B matrix
        E % E matrix
        C % C matrix
        G % G matrix
        n % state dimension
        m % input dimension
        l % disturbance dimension
        k % output dimension
        r % noise dimension
    end
    
    methods
        function obj = LTISys(A,B,E,C,G)
            
            obj.A = A;
            obj.B = B;
            obj.E = E;
            obj.C = C;
            obj.G = G;
            
            obj.n = size(A,1);
            obj.m = size(B,2);
            obj.k = size(C,1);
            obj.l = size(E,2);
            obj.r = size(G,2);
            obj.assertValidSystem();
            
        end
        
        function [] = assertValidSystem(obj)
            assert(all(size(obj.A) == [obj.n,obj.n]));
            assert(all(size(obj.B) == [obj.n,obj.m]));
            assert(all(size(obj.C) == [obj.k,obj.n]));
            assert(all(size(obj.E) == [obj.n,obj.l]));
            assert(all(size(obj.G) == [obj.k,obj.r]));
        end
        
        function [Ky,uc] = computeOutputFeedbacks(obj,horizon,Kyt,uct)
            toe_mat = [zeros(obj.k,obj.m * (horizon-1));
                       LTISys.computeToeplitzMat(obj.A,obj.B,obj.C,horizon-1)];
            toe_mat = [toe_mat, zeros(size(toe_mat,1),obj.m)];
            Ky = linsolve(eye(size(Kyt,1))+Kyt*toe_mat, Kyt);
            uc = uct;
        end
        function [Kyt,uct] = computeModifiedOutputFeedbacks(obj,horizon,Ky,uc)
            toe_mat = [zeros(obj.k,obj.m * (horizon-1));
                       LTISys.computeToeplitzMat(obj.A,obj.B,obj.C,horizon-1)];
            toe_mat = [toe_mat, zeros(size(toe_mat,1),obj.m)];
            Kyt = linsolve(eye(size(Ky,2))-Ky'*toe_mat', Ky')';
            uct = uc;
        end
        
        function [sysmap, x_tmat, w_tmat] = computeSystemMapModified(obj, horizon, Kyt, uct)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obs_tmat = LTISys.computeToeplitzMat(obj.A,eye(obj.n),obj.C,horizon);
            Uxw = Kyt * [obs_tmat,zeros(horizon * obj.k,obj.n)];
            Ux = Uxw(:,1:obj.n);
            
            Uw = Uxw(:,obj.n+1:end) * kron(eye(horizon),obj.E);
            Uv = Kyt * kron(eye(horizon),obj.G);
            u = uct;
            
            xw_tmat = LTISys.computeToeplitzMat(obj.A,eye(obj.n),eye(obj.n),horizon+1);
            x_tmat = xw_tmat(:,1:obj.n);
            w_tmat = xw_tmat(:,obj.n+1:end) * kron(eye(horizon),obj.E);
            u_tmat = [zeros(obj.n,horizon * obj.m); LTISys.computeToeplitzMat(obj.A,obj.B,eye(obj.n),horizon)];
            Sx = x_tmat + u_tmat * Ux;
            Sw = w_tmat + u_tmat * Uw;
            Sv = u_tmat * Uv;
            s = u_tmat * u;
            
            sysmap.Sx = Sx;
            sysmap.Sw = Sw;
            sysmap.Sv = Sv;
            sysmap.s = s;
            sysmap.Ux = Ux;
            sysmap.Uw = Uw;
            sysmap.Uv = Uv;
            sysmap.u = u;
            sysmap.horizon = horizon;
        end
        
        function [sysmap, x_tmat, w_tmat] = computeSystemMap(obj,horizon, Ky, uc)
            [Kyt,uct] = obj.computeModifiedOutputFeedbacks(horizon,Ky,uc);
            [sysmap, x_tmat, w_tmat] = computeSystemMapModified(obj, horizon, Kyt, uct);
        end
        
        function [x,u,w,v] = simulateClosedLoopRealization(obj,sysmap,x0,wcol,vcol)
            xcol = sysmap.Sx * x0 + sysmap.Sw * wcol + sysmap.Sv * vcol + sysmap.s;
            ucol = sysmap.Ux * x0 + sysmap.Uw * wcol + sysmap.Uv * vcol + sysmap.u;
            
            x = reshape(xcol,[obj.n,sysmap.horizon+1]);
            u = reshape(ucol,[obj.m,sysmap.horizon]);
            w = reshape(wcol,[obj.l,sysmap.horizon]);
            v = reshape(vcol,[obj.r,sysmap.horizon]);
        end
        
    end
    
    methods(Static)
        
        function Tmat = computeToeplitzMat(A,B,C,N)
            if N == 0
                Tmat = zeros(0);
                return
            end
            blocks = cell(N+1,1);
            blocks{1} = zeros(size(C,1),size(B,2));
            for i = 1:N
                blocks{i+1} = C * A^(i-1) * B;
            end
            Tmat = cell2mat(blocks(tril(toeplitz(1:N))+1)); 
        end
    end
end

