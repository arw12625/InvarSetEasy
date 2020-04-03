function [q, converged, infeasible] = insideOutRecurrence(sys,cf,uncertainty,admissible, P, q0, c, max_iter, epsilon, yalmipOptions)
%COMPUTEAFFINECONTROLLER Summary of this function goes here
%   Detailed explanation goes here

q = q0;
qvar = sdpvar(size(q,1),1,'full');

A_un = [uncertainty.A; uncertainty.Ae; -uncertainty.Ae;
        [P, zeros(size(P,1),size(uncertainty.A,2) - size(P,2))]];
A_ad = [admissible.A; admissible.Ae; -admissible.Ae;
        [zeros(size(P,1),size(admissible.A,2) - size(P,2)), P]];
b_ad = [admissible.b; admissible.be; -admissible.be;
        qvar - epsilon * ones(size(q,1),1)];

sysmap = cf.computeSystemMap();

F = [sysmap.Ux,sysmap.Uw,sysmap.Uv
     sysmap.Sx,sysmap.Sw,sysmap.Sv;];
f = [sysmap.u;
     sysmap.s];

T = sdpvar(size(A_ad,1), size(A_un,1),'full'); 

converged = 0;
infeasible = 0;
for iter = 1:max_iter
    iter
    b_un = [uncertainty.b; uncertainty.be; -uncertainty.be;
            q];

        
    constraints = [cf.constraints;
                   T >= 0;
                   T * A_un == A_ad * F;
                   T * b_un <= b_ad - A_ad * f;
                   qvar >= epsilon * ones(size(q,1),1);
                   qvar >= q * 0.9];

    diagnostics = optimize(constraints, c'*qvar, yalmipOptions);
    
    if diagnostics.problem ~= 0
        yalmiperror(diagnostics.problem)
        infeasible = 1;
        break
    end
    q_old = q;
    q = value(qvar);
    %c'*q
    q'
    if all(q <= q_old + epsilon * ones(size(q,1),1))
        converged = 1;
        break
    end
end
end

