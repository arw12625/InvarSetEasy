function [p_ineq] = eq2ineq(p_eq)
%EQ2INEQ Transforms a polyhedron with equality constraints into one with
%only inequality constraints.

A_ineq = [p_eq.A;
          p_eq.Ae;
          -p_eq.Ae;];
b_ineq = [p_eq.b;
          p_eq.be;
          -p_eq.be;];
p_ineq = Polyhedron(A_ineq, b_ineq);

end