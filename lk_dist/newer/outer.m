%%

scale = 1;
N = 5;

%% lin method

s = scale * Omega_lk;
arr = [];
for i = 1:N
    s = polyLinPre(plsys, s);
    arr = [arr, s];
end
polyunion = PolyUnion(arr);
inner_invar = polyunion.convexHull();
volume(inner_invar)

%% proj method
s = X_lk;

steps = 10;

start_time = cputime;
for i = 1:steps
    s = s & polyLinPre(plsys, s);
    [i, size(s.A,1)]
end

end_time = cputime - start_time
outer_invar = s;
%volume(outer_invar)

