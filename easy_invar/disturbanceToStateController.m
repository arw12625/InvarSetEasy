function [K,uc] = disturbanceToStateController(plsys, N, Kx,Kw,u0)
%disturbanceToStateController Compute a state feedback controller
%corresponding to the given disturbance feedback controller
%   plsys - the system to be controlled
%   N - number of steps the controller determines input for
%   Kx,Kw,u0 - the disturbance feedback controller provided
%       u = Kx x_0 + Kw w + u0
%       where x_0 is the initial state and w is the disturbance sequence
%       note that u,w are in reverse order, the last input is first
%
%   K,uc - the state feedback controller found
%       u = K x + uc
%       note here x is the state trajectory, and both u,x are specified in
%       forward order unlike the disturbance controller

n = plsys.n;
m = plsys.m;
l = plsys.l;

Einv = pinv(plsys.E);

Kx = mat2cell(Kx, repmat(m,1,N),n);
Kw = mat2cell(Kw, repmat(m,1,N),repmat(l,1,N));
u0 = mat2cell(u0, repmat(m,1,N),1);

K = cell(N,1);
uc = cell(N,1);

K{1} = [Kx{N},zeros(m, n*(N-1))];
uc{1} = u0{N};

for i = 2:N
    M = zeros(m, N * n);
    M(:,1:n) = Kx{N - i + 1};
    uc{i} = u0{N+1-i};
    for j = 2:i
        M(:,(j-1)*n+(1:n)) = M(:,(j-1)*n+(1:n)) + Kw{N+1-i,N+1-j+1} * Einv;
        M(:,(j-2)*n+(1:n)) = M(:,(j-2)*n+(1:n)) - Kw{N+1-i,N+1-j+1} * Einv * plsys.A;
        M = M - Kw{N+1-i,N+1-j+1} * Einv * plsys.B * K{j-1};
        uc{i} = uc{i} - Kw{N+1-i,N+1-j} * Einv * plsys.B * u0{N+1-j};
    end
    K{i} = M;
end

K = cell2mat(K);
uc = cell2mat(uc);

end

