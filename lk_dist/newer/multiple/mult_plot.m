 
s = Omega_lk;

for i = 1:N
    s = polyLinPre(plsys, s);
end

%{
co = [222,235,247] / 255;
cs = [158,202,225] / 255;
cx = [49,130,189] / 255;
%}

co = [0.9, 0.3, 0.3];
cs = [0.3, 0.9, 0.3];
cx = [0.3, 0.3, 0.9];

figure
plot(slice(X_lk,[2,3,4,6,7,8]), 'color', cx, slice(s, [2,3,4,6,7,8]), 'color', cs, slice(Omega_lk, [2,3,4,6,7,8]), 'color', co);
legend('X_s', 'Pre^N(\Omega)', '\Omega');
legend('location', 'southeast');
title('Lane Keeping Polytope Slice corresponding to v_i = \Psi_i = r_i = 0 for i=1,2')