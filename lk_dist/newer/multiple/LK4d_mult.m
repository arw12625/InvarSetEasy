%Constructs LK 4D Matrices 

  %Params

  M = 1800;
  a = 1.2;
  b = 1.65;
  Caf = 140000;
  Car = 120000;
  Iz = 3270;
  u = 10;

  %Dynamics
  % x(k+1) = Ax(k) + Bu(k) + Ew(k) + f
  % x = [y; v; psi; r]
  % u = df;
  % w = rd;
  
  A_lk = [0, 1, u, 0; 
          0, -(Caf+Car)/M/u, 0, ((b*Car-a*Caf)/M/u - u); 
          0, 0, 0, 1;
          0, (b*Car-a*Caf)/Iz/u,  0, -(a^2 * Caf + b^2 * Car)/Iz/u];
  B_lk = [0; Caf/M; 0; a*Caf/Iz];
  E_lk = [0; 0; -1; 0];
  
  csys = ss(A_lk, [B_lk, E_lk], eye(size(A_lk,2)), zeros(size(A_lk, 2), size([B_lk,E_lk],2)));
  dsys = c2d(csys, 0.2);
  A_lk = dsys.A;
  B_lk = dsys.B(:,1:size(B_lk, 2));
  E_lk = dsys.B(:,size(B_lk, 2) + (1:size(E_lk, 2)));
  f_lk = zeros(size(A_lk,1),1);
  
  d = 2; %number of vehicles
  
  A_lk = kron(eye(d), A_lk);
  B_lk = kron(eye(d), B_lk);
  E_lk = repmat(E_lk, d, 1);
  f_lk = repmat(f_lk, d, 1);
  
  %Safe Set
  
    %Lat vel    +/-     1.0   m/s
    %yaw rate   +/-     0.78  rad/s
    %Offset     +/-     0.7   m
   
   yH = 2;
   yL = -2;
   vH = 1;
   vL = -1;
   rH = 0.78;
   rL = -0.78;
   %{
   dfH = 0.7854;
   dfL = -0.7854;
   rdH = 0.5;
   rdL = -0.5;
   %}
   
   dfH = 1;
   dfL = -1;
   rdH = 0.3;
   rdL = -0.3;
   
   Axh_lk =  [1, 0, 0 ,0;
         -1, 0, 0, 0;
         0, 1, 0 ,0;
         0, -1, 0, 0;
         0, 0, 0 ,1;
         0, 0, 0, -1];
   bxh_lk = [yH;-yL;vH;-vL;rH;-rL];
   
   %sep_dist = 0.25;
   sep_dist = 1.5;
   
   X_lk = Polyhedron(kron(eye(d), Axh_lk), repmat(bxh_lk, d, 1));
   
   sep_poly = Polyhedron([1, 0, 0, 0, -1, 0, 0, 0], -sep_dist);
   
   X_lk = X_lk & sep_poly;
   
  %Constraints
  
    %driving angle +/-45pi/180
    
  U_lk = Polyhedron( 'H', ...
        [1, 0, dfH;
         -1, 0, -dfL
         0, 1, dfH;
         0, -1, -dfL]);
     
     
  W_lk = Polyhedron( 'H', ...
        [1, rdH;
         -1, -rdL]);
  
     %{
  Omega_lk = Polyhedron( 'H', ...
        [1, 0, 0 ,0, yH / 3;
         -1, 0, 0, 0, -yL / 3;
         0, 1, 0 ,0, vH / 2.5;
         0, -1, 0, 0, -vL / 2.5;
         0, 0, 1 ,0, rH / 5;
         0, 0, -1, 0, -rL / 5;
         0, 0, 0 ,1, rH / 3;
         0, 0, 0, -1, -rL / 3;]);
         %}
         Omega_lk = Polyhedron( 'H', ...
        [1, 0, 0 ,0, yH / 5;
         -1, 0, 0, 0, -yL / 5;
         0, 1, 0 ,0, vH / 5;
         0, -1, 0, 0, -vL / 5;
         0, 0, 1 ,0, rH / 5;
         0, 0, -1, 0, -rL / 5;
         0, 0, 0 ,1, rH / 5;
         0, 0, 0, -1, -rL / 5;]);
    
     offset = sep_dist / 2 + yH / 5 + yH / 10;
   Omega_lk = ([-offset;0;0;0] + Omega_lk) * ([offset;0;0;0] + Omega_lk);
   %Omega_lk = Omega_lk & 3 * sep_dist;

plsys = PolyLinSys(A_lk, X_lk, B_lk, U_lk, E_lk, W_lk, f_lk);
pslsys = PolySwitchLinSys({A_lk},X_lk,{B_lk},U_lk,{E_lk}, W_lk, {f_lk});