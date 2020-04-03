clear all
% Dynamics setting for TCLs

N = 10000 ;     % Number of TCLs
th_a = 32 ;     % Ambient temperature

%% TCL's property constants

c = 3 ;         
r = 2 ;
P = 5.6 ;
eta = 2.5 ;
th_r = 22.5 ;
del_th = 1.2 ;

tcl_1 = tcl(c,r,P,eta,th_r,del_th) ;

%% Abstraction parameters 

tau = 0.2 ;         % Time discretization 
%eta = 0.05 ;         % State(temperature) discretization
%eps = 0.05 ;         % Bisimilarity ( should be bigger than eta/2)
eta = 0.2;
eps = 0.2;

%% Generate system matrices 

[A,B,K,lowerBound,upperBound,nd_upp,nd_low] = systemPetter(tcl_1,th_a,N,tau,eta) ;     % Generate abstraction model

%% System constraints as polyhedra
%{
% N = 100;                  % number of loads
upperBound = 9000;          % upper bound on the number of loads that are on
lowerBound = 5000;          % upper bound on the number of loads that are on

% This polyhedron models that constraints for the upper and lower bound on
% the number of active loads, the nonnegativity requirements, and the total
% number of loads.
powerConstraints = Polyhedron('H', ...
   [ones(1,K), zeros(1,K), N - lowerBound;
    -eye(K), zeros(K), zeros(K,1);
    zeros(1,K), ones(1,K), upperBound;
    zeros(K), -eye(K), zeros(K,1)], 'He', ...
    [ones(1,2*K), N] );

% This polyhedron models the constraints that the number of loads to switch
% from m1 to m2 cannot be more than the number of loads at m1, and that the
% number of loads switched must be nonnegative.
jointInputConstraints = Polyhedron('H', ...
    [-eye(K), zeros(K), eye(K), zeros(K), zeros(K,1);
     zeros(K), zeros(K), -eye(K), zeros(K), zeros(K,1);
     zeros(K), -eye(K), zeros(K), eye(K), zeros(K,1);
     zeros(K), zeros(K), zeros(K), -eye(K), zeros(K,1)]);
 
% By the condition of unsafe set of abstraction in paper () 
% we add that the number of TCLs in the nodes next to the boundary should
% be zero.

fullInputSpace = Polyhedron('H', [zeros(1,2*K), 1]);
                        

% This polyhedron is the combination of the state and input constraints
XU = (powerConstraints * fullInputSpace) & jointInputConstraints;
Xterm = powerConstraints;

% For now we ignore disturbance
W = Polyhedron();
%}