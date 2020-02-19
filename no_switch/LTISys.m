classdef LTISys < matlab.mixin.Copyable
    
    properties (SetAccess = immutable)
        A % A matrix
        B % B matrix
        C % C matrix
        n % state dimension
        m % input dimension
        k % output dimension
    end
    
    methods
        function obj = LTISys(A,B,C)
            
            obj.A = A;
            obj.B = B;
            obj.C = C;
            
            obj.n = size(A,1);
            obj.m = size(B,2);
            obj.k = size(C,1);
            obj.assertValidSystem();
            
        end
        
        function [] = assertValidSystem(obj)
            assert(all(size(obj.A) == [obj.n,obj.n]));
            assert(all(size(obj.B) == [obj.n,obj.m]));
            assert(all(size(obj.C) == [obj.k,obj.n]));
        end
        
    end
    
    methods(Static)
        
        function Tmat = computeToeplitzMat(A,B,C,N)
            blocks = cell(N+1,1);
            blocks{1} = zeros(size(C,1),size(B,2));
            for i = 1:N
                blocks{i+1} = C * A^(i-1) * B;
            end
            Tmat = cell2mat(blocks(tril(toeplitz(1:N))+1)); 
        end
    end
end

