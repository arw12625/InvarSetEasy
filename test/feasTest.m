%Test control invariance for a given system with state constraints

Omega = 1 / 0.81 * Polyhedron.unitBox(2); %Initial Set
X = 5 * Polyhedron.unitBox(2); %State Set
U = 0.6 * Polyhedron.unitBox(2); %Input Set

theta = pi / 4;
A = 1.4 * [cos(theta), sin(theta); -sin(theta), cos(theta)];
B = eye(2);

N = 1; %Number of backwards steps to use


Omega = Polyhedron.unitBox(n);
X = 10 * Polyhedron.unitBox(n);
U = Polyhedron.unitBox(n);

A = 2 * [1,0;0,1];
B = [1,1;0,1];

testStateControlInvariance(Omega, X, U, N, A, B)

%{

%N = 2;
s = Omega;
clf;
hold on;

for N = 1:5
    plot(s); 
    s = polyPre(A,B,X,U,s);
    testStateControlInvariance(Omega, X, U, N, A, B)
end
alpha(0.1);

%}


