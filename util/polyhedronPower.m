function [pn] = polyhedronPower(poly, n, irredundant)
%POLYHEDRONPOWER Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    irredundant = false;
end

A = kron(eye(n), poly.A);
b = repmat(poly.b, n,1);
Ae = kron(eye(n), poly.Ae);
be = repmat(poly.be, n,1);

pn = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be, 'irredundantHRep', irredundant);

end

