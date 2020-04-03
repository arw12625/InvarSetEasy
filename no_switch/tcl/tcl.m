classdef tcl
    %TCL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c ;
        r ;
        P ;
        eta ;
        th_r ;
        del_th ;
        a ;
        b ;
        lb ;
        ub ;
        
    end
    
    methods
        
        function obj = tcl(c,r,P,eta,th_r,del_th)
            %TCL Construct an instance of this class
            %   Detailed explanation goes here
            obj.c = c ;
            obj.r = r ;
            obj.P = P ;
            obj.eta = eta ;
            obj.a = 1 / (c*r) ;
            obj.b = eta / c ;
            obj.th_r = th_r ;
            obj.del_th = del_th;
            obj.lb = th_r - del_th ;
            obj.ub = th_r + del_th ;
            
        end
        
    end
    
end

