function [A,B,K,max_low,min_upp,uppnodes,lownodes,coff,con] = systemPetter(tcl,th_a,N,tau,eta)

    %% Parameters of TCLs
    
    a = tcl.a ;     
    b = tcl.b ;
    Pm = tcl.P ;
    
    lb = tcl.lb ;
    ub = tcl.ub ;

    max_low = N*a*(th_a - lb)/ (b*Pm) ;
    min_upp = N*(1 - a*(th_a - ub)/(b*Pm)) ;

    %% Number of the nodes

    K = ceil( (ub - lb) / eta );

    % Compute the lower, upper, and mid temperature of each state

    n_lb = linspace(lb, ub-eta, K) ;
    n_ub = linspace(lb+eta, lb+K*eta, K) ;

    if n_ub(end) > ub
        n_ub(end) = ub; 
    end

    n_mid = mean([n_lb ; n_ub]) ;
    
    % Generate System Matrices

    A = zeros(2*K,2*K) ;          % State transition matrix
    B = zeros(2*K,2*K) ;          % Input transition matrix
    
    uppnodes = [] ;
    lownodes = [] ;

    for i = 1 : K

        n1 = temp2node(n_mid(i),lb,eta) ;
        nx_temp = nexttemp(n_mid(i),0,tcl,tau,th_a) ;
        n20 = temp2node(nx_temp,lb,eta) ;
        
        if n1 == n20 
            
            fprintf("Repetitive node exists!") ;
            assert(n1 ~= n20) ;
            
        end

        if n20 <= K

            A(n20,n1) = 1 ;
            B(n20,n1) = -1 ;
            B(n20,n1+K) = 1 ;
            
        elseif n20 > K
            
            uppnodes = [uppnodes ; n1] ;
            
        end
        
        if i == 1
            
            coff = n20 ;
            
        end

        n1 = temp2node(n_mid(i),lb,eta) ;
        nx_temp = nexttemp(n_mid(i),1,tcl,tau,th_a) ;
        n20 = temp2node(nx_temp,lb,eta) ;

        if n20 > 0

            A(n20+K,n1+K) = 1 ;
            B(n20+K,n1) = 1 ;
            B(n20+K,n1+K) = -1 ;
            
        elseif n20 <= 0
            
            lownodes = [lownodes ; n1] ;

        end
        
        if i == K
            
            con = n20 ;
            
        end

    end

end


function node = temp2node(temp,lb,eta)

    node = floor((temp - lb) / eta) + 1 ; 

end

function n_temp = nexttemp(temp,mode,tcl,tau,th_a)

    a = tcl.a ;
    b = tcl.b ;
    Pm = tcl.P ;

    if mode == false
        
        n_temp = temp*exp(-a*tau) + th_a * (1 - exp(-a*tau)) ;
        
    elseif mode == true
        
        n_temp = temp*exp(-a*tau) + (th_a - b*Pm/a) * (1 - exp(-a*tau)) ;
        
    end
    
end
