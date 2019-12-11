
sbv = @(i,d) [zeros(i-1,1);1;zeros(d-i,1)];
T = 5;

%{
d = 5;
adj_pairs = [1,2;
       2,3;
       3,1;
       1,4;
       4,5;
       5,5];
%}

d = 2;
adj_pairs = [1,1;2,2];

num_pairs = size(adj_pairs,1);
XU_vecs = zeros(num_pairs, 2*d);
for i = 1:num_pairs
    XU_vecs(i,:) = [sbv(adj_pairs(i,1),d)',sbv(adj_pairs(i,2),d)'];
end
XU = Polyhedron('V',XU_vecs);

A = zeros(d);
B = eye(d);
E = zeros(d,0);
f = zeros(d,1);

W = Polyhedron([],[]);

Xterm = Polyhedron(zeros(1,d),1);
Xterm = eq2ineq(Xterm);

sys = LTISys(T,A,B,E,f,XU,Xterm,W);

init = Polyhedron('V',sbv(1,d)');
init = eq2ineq(init);
init = Xterm;

[isRecurrent, affineController] = testAffineRecurrence(sys,init);