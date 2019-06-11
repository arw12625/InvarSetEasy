classdef PolyLinSys
    %PolyLinSys represents a linear system with polytopic regions
    %The dynamics of the system are given by 
        %x(k+1) = Ax(k) + Bu(k) + Ew(k) + f
        % x \in X, u \in U, w \in W
    %X,U,W are all MPT polytopes
    
    properties
        A
        B
        E
        f
        X
        U
        W
        n % state dimension
        m % input dimension
        l % disturbance dimension
    end
    
    methods
        function obj = PolyLinSys(A,X,B,U,E,W,f)
            %PolyLinSys Construct an instance of this class
            if nargin == 7 || nargin == 6
                obj.A = A;
                obj.B = B;
                obj.E = E;
                obj.X = X;
                obj.U = U;
                obj.W = W;
                obj.n = size(A,1);
                obj.m = size(B,2);
                obj.l = size(E,2);
                if nargin == 6
                    obj.f = zeros(obj.n,1);
                else 
                    obj.f = f;
                end
            elseif nargin == 4
                obj.A = A;
                obj.B = B;
                obj.X = X;
                obj.U = U;
                obj.n = size(A,1);
                obj.m = size(B,2);
                obj.l = 0;
                obj.E = zeros(size(obj.n, obj.l));
                obj.W = Polyhedron;
                obj.f = zeros(obj.n,1);
            else
                error('Wrong number of input arguments')
            end
        end
        
        function obj = createEmptySys(n,m,l)
            %createEmptySys initialize a system with the given dimensions
                %with zero dynamics and empty regions
            obj = PolyLinSys( ...
                zeros(n), ...
                Polyhedron('H',zeros(1,n+1)), ...
                zeros(n,m), ...
                Polyhedron('H',zeros(1,m+1)), ...
                zeros(n,l), ...
                Polyhedron('H',zeros(1,l+1)), ...
                zeros(n,1));
        end
        
        function dist = hasDisturbance(obj)
            %hasDisturbance determine if the system has a nontrivial disturbance
            dist = (obj.E ~= 0) && ~(isEmptySet(obj.W));
        end
    end
end

