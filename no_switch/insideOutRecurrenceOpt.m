function [q, converged, infeasible] = insideOutRecurrenceOpt(sys,cf,uncertainty,admissible, P, q0, c, max_iter, epsilon, yalmipOptions)
%COMPUTEAFFINECONTROLLER Summary of this function goes here
%   Detailed explanation goes here

q_old = sdpvar(size(q0,1),1,'full');
q_new = sdpvar(size(q0,1),1,'full');

A_un = [uncertainty.A; uncertainty.Ae; -uncertainty.Ae;
        [P, zeros(size(P,1),size(uncertainty.A,2) - size(P,2))]];

A_ad = [admissible.A; admissible.Ae; -admissible.Ae;
        [zeros(size(P,1),size(admissible.A,2) - size(P,2)), P]];

b_un = [uncertainty.b; uncertainty.be; -uncertainty.be;
        q_old];

b_ad = [admissible.b; admissible.be; -admissible.be;
        q_new - epsilon * ones(size(q0,1),1)];

T = sdpvar(size(A_ad,1), size(A_un,1),'full');

sysmap = cf.computeSystemMap();

F = [sysmap.Ux,sysmap.Uw,sysmap.Uv
     sysmap.Sx,sysmap.Sw,sysmap.Sv;];
f = [sysmap.u;
     sysmap.s];

alpha = sdpvar(1);
base_alpha = 0.9;
alphas = [base_alpha .^ ((floor(max_iter/2)-1):-1:0), ones(1,ceil(max_iter/2))];
alphas = max(alphas,epsilon);

gamma = sdpvar(1);

constraints = [cf.constraints;
               T >= 0;
               T * A_un == A_ad * F;
               T * b_un <= b_ad - A_ad * f;
               q_new >= epsilon * ones(size(q0,1),1);
               q_new >= q_old * alpha;
               q_new <= gamma * ones(size(q0,1),1)];

objective = c' * q_new + gamma;

P_optimizer = optimizer(constraints, objective, yalmipOptions, {q_old, alpha}, q_new);

q = q0;
obj_val = c'*q;
converged = 0;
infeasible = 0;
for iter = 1:max_iter
    iter

    [qt, yalmip_num] = P_optimizer(q,alphas(iter));
    %alphas(iter)
    %c'*qt
    %qt'
    if yalmip_num ~= 0
        yalmiperror(yalmip_num)
        infeasible = 1;
        break
    end
    if all(qt <= q + epsilon * ones(size(q,1),1))
        q = qt;
        converged = 1;
        break
    end
    obj_val = c'*qt
    q = qt;
    max(q')
end
end

